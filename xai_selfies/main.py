import selfies as sf
from collections import defaultdict
import numpy as np
import pandas as pd
import re
from tqdm import tqdm
from rdkit.Chem import MolFromSmiles
from rdkit import Chem, DataStructs


NON_ATOM_CHARACTERS = ['</s>', '#', '%', ')', '(', '+', '-', '1', '0', '3', '2', '5', '4', '7', '6', '9', '8', ':',
                       '=', '@', '[', ']', '<s>']


def smiles_parser(sml):
    """Parse the input smiles to get list of tokens (char)
    :param sml: SMILES string
    :return list of tokens that are part of the SMILES string, indices of tokens that correspond to atoms, list
    containing only the tokens at the positions that are atoms
    """
    REGEX_SML = r'Cl|Br|[#%\)\(\+\-1032547698:=@CBFIHONPS\[\]cionps]'
    atom_chars = ["Cl", "Br", "C", "B", "F", "I", "H", "O", "N", "P", "S", "c", "i", "o", "n", "p", "s"]
    char_list = re.findall(REGEX_SML, sml)
    atom_char_list = [i for i in range(len(char_list)) if char_list[i] in atom_chars]
    true_char_list = [char_list[i] for i in atom_char_list]
    return char_list, atom_char_list, true_char_list


def get_correct_order(input_smiles):
    """
    The SMILES given by the user might not be RDKit canonical SMILES. The trick is to re-order the attribution
    weights to the RDKit canonical atom ordering to make sure all things are comparable.
    :param input_smiles: input SMILES from whichever source
    :return: numpy array with the indices to reorder the attributions to follow the atom ordering from RDKit
    canonicalization algorithm
    """
    mol = Chem.MolFromSmiles(input_smiles)
    Chem.MolToSmiles(mol, canonical=True)
    atom_order = mol.GetProp("_smilesAtomOutputOrder")
    atom_order = atom_order[1:-2].split(',')
    return np.array([int(ao) for ao in atom_order])


def get_mutated_selfies(selfie, position, similarity=0.0):
    """
    Given a SELFIES string and a position that we want to "mutate", returns all possible
    alternative SELFIES that correspond to valid molecules (optionally within a similarity threshold of the query)
    :param selfie: SELFIES string for the query molecule
    :param position: which index in the string must be mutated
    :param similarity: minimum Tanimoto similarity needed to count a mutated molecule as a parent of
    the query molecule. Default 0.0, in which case no filtering based on similarity will be performed. Valid
    similarities are floats between 0 and 1.0, 0.0 allowing all possible mutants and 1.0 effectively filtering out all
    mutated molecules. A recommended value could be around 0.5. Note: any value different from 0.0 will decrease the
    speed of the algorithm significantly.
    :return: a list of mutated selfies and a list of the corresponding SMILES strings
    """
    mutated_selfies = []
    mutated_smiles = []
    mol = Chem.MolFromSmiles(sf.decoder(selfie))
    chars_selfie = list(sf.split_selfies(selfie))
    alphabet = list(sf.get_semantic_robust_alphabet())  # 69 SELFIE tokens
    for character in alphabet:
        selfie_mutated_chars = chars_selfie[:position] + [character] + chars_selfie[position + 1:]
        selfie_mutated = "".join(x for x in selfie_mutated_chars)
        # Check for molecular validity
        try:
            smiles = sf.decoder(selfie_mutated)
            rdkit_mol = MolFromSmiles(smiles, sanitize=True)
            if rdkit_mol is not None:
                mutated_selfies.append(selfie_mutated)
                mutated_smiles.append(smiles)
        except sf.DecoderError:
            print('Invalid molecule!')
            pass
        if similarity > 0.:
            # otherwise no need to bother with filtering
            assert 1.0 >= similarity >= 0.0
            mutated_smiles, close_indices = filter_candidates(mol, mutated_smiles, similarity)
            mutated_selfies = list(np.array(mutated_selfies)[close_indices])

    return mutated_selfies, mutated_smiles


def filter_candidates(mol, mutated_smiles, similarity=0.5):
    """
    Keep from the mutated molecules only those within a similarity threshold.
    :param mol: RDKit molecule to compare our mutated molecules with
    :param mutated_smiles: list of SMILES strings to filter for dissimilar molecules
    :param similarity: minimum fingerprint similarity the mutated molecules are allowed to have with respect to the
    query mol
    :return: list of SMILES that pass the similarity threshold, list of indices to the original list
    """
    close_smiles = []
    close_indices = []
    mols = [Chem.MolFromSmiles(smi) for smi in mutated_smiles]
    fps = [Chem.RDKFingerprint(x) for x in mols]
    original_fps = Chem.RDKFingerprint(mol)
    for i in range(len(mutated_smiles)):
        sim = DataStructs.FingerprintSimilarity(original_fps, fps[i])
        if sim >= similarity:
            close_smiles.append(mutated_smiles[i])
            close_indices.append(i)
    return close_smiles, close_indices


def score_from_selfies_to_smiles(selfie, scores):
    """
    Given a SELFIES string and the computed attribution scores at the selfie token level, convert this to a SMILES
    attribution scores (at the SMILES character level) using the Attribute functionality in sf.
    :param selfie: SELFIES string corresponding to one molecule
    :param scores: list of attribution scores (len(selfie_tokens))
    :return: list of attribution scores (len(smiles_tokens))
    """
    smiles, attributions = sf.decoder(selfie, attribute=True)
    smiles_tokens, atom_char_list, true_char_list = smiles_parser(smiles)
    skipped_chars = ['(', ')', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ']', '[', '@', '%']
    smiles_scores = defaultdict(float)
    smiles_scores_list = []
    for smiles_token in attributions:
        if smiles_token.attribution is None:  # can happen, not sure why
            smiles_scores[str(smiles_token.index) + '_' + smiles_token.token] = 0.0
        else:
            for att in smiles_token.attribution:
                if '[' in smiles_token.token:  # cases like [NH1]:
                    for char in smiles_token.token[1:-1]:
                        if char not in skipped_chars:
                            smiles_scores[str(smiles_token.index) + '_' + char] += scores[att.index]
                else:
                    smiles_scores[str(smiles_token.index) + '_' + smiles_token.token] += scores[att.index]
    # add zeros for missing smiles tokens in the translation
    i = 0
    smiles_weights = list(smiles_scores.values())
    for smiles_char in smiles_tokens:
        if smiles_char in skipped_chars:
            smiles_scores_list.append(0.0)
        else:
            smiles_scores_list.append(smiles_weights[i])
            i += 1
    assert len(smiles_tokens) == len(smiles_scores_list)

    # Keep only weights for atoms
    atom_characters_indices = [i for i, character in enumerate(smiles_tokens) if character not in NON_ATOM_CHARACTERS]
    atom_weights = np.array(smiles_scores_list)[atom_characters_indices]
    return selfie, scores, smiles, smiles_tokens, smiles_scores_list, atom_weights


def get_all_mutations(smiles, ids=None, similarity=0.):
    """

    :param smiles: list of SMILES string for which we want to obtain all possible SELFIES mutations
    :param ids: if given, list of unique identifiers for the molecules
    :param similarity: minimum fingerprint similarity the mutated molecules are allowed to have with respect to the
    query mol
    :return: dataframe with original smiles, selfies, and mutations
    """
    original_smiles = []
    mutated_smiles = []
    original_selfies = []
    mutated_selfies = []
    mutation_position = []
    molids = []
    if isinstance(smiles, list):
        for i in tqdm(range(len(smiles))):
            sml = smiles[i]
            print(f"Mutating {sml}")
            orig_selfie = sf.encoder(sml)
            chars = list(sf.split_selfies(orig_selfie))
            for j, sel_char in enumerate(chars):
                mutated_self, mutated_smls = get_mutated_selfies(orig_selfie, position=j, similarity=similarity)
                original_smiles += [sml for _ in mutated_self]
                original_selfies += [orig_selfie for _ in mutated_self]
                mutated_smiles += mutated_smls
                mutated_selfies += mutated_self
                mutation_position += [j for _ in mutated_self]
                if ids is None:
                    molids += [i for _ in mutated_self]
                else:
                    molids += [ids[i] for _ in mutated_self]
        return pd.DataFrame.from_dict({'molid': molids, 'smiles': original_smiles, 'selfies': original_selfies,
                                       'mutated_smiles': mutated_smiles, 'mutated_selfies': mutated_selfies,
                                       'mutation_pos': mutation_position})
    else:
        print('Wrong input format, expected list of SMILES')
        return


def get_predictions_for_mutants_and_original(mutant_df, model, embedder):
    """

    :param mutant_df: dataframe containing original smiles, mutated smiles and positions (obtained from
    get_all_mutations)
    :param model: python object with a predict() method returning continuous-valued predictions (regression
    model or target class probability) given an adequately featurized input
    :param embedder: python object with an encode() method able to take as input a list of SMILES strings and
    returning the expected input features for trained_model
    :return: dataframe with additional column with predictions from the model, array of original predictions
    """
    unique_orig = mutant_df.drop_duplicates(subset=['molid', 'smiles'])['smiles'].tolist()
    mutated = mutant_df['mutated_smiles'].tolist()
    all_smis = unique_orig + mutated
    X = np.array(embedder.encode(all_smis))
    preds = model.predict(X)
    mutant_df['predictions'] = preds[len(unique_orig):]
    return mutant_df, preds[:len(unique_orig)]


def get_attributions_df(smiles_list, model, embedder, ids=None, similarity=0.):
    """
    Given a list of molecules represented as SMILES string, a trained model to explain, a feature extractor, and a
    similarity threshold, compute the SELFIES atom attribution method and return results as a pandas dataframe.
    :param smiles_list: list of input molecules for which we want to obtain predictions and explanations as atom
    attributions
    :param model: trained model to explain, python object with a predict() method returning continuous-valued
    predictions (regression model or target class probability) given an adequately featurized input.
    :param embedder: python object with an encode() method able to take as input a list of SMILES strings and
    returning the expected input features for the model
    :param ids: optional: unique identifier for the molecules
    :param similarity: minimum fingerprint similarity the mutated molecules are allowed to have with respect to the
    query molecule
    :return: pandas dataframe with one row per input smiles. Returns all necessary properties including predictions and
    atom attributions.
    """

    if not hasattr(model, 'predict'):
        raise AttributeError('The model must implement a predict() method.')
    if not hasattr(embedder, 'encode'):
        return AttributeError('The featurizer must implement an encode() method.')

    # 1. Get all mutants
    mutant_df = get_all_mutations(smiles_list, ids=ids, similarity=similarity)

    # 2. Get all predictions
    mutant_df, original_preds = get_predictions_for_mutants_and_original(mutant_df, model=model, embedder=embedder)
    print('Got all mutants for the query molecules')
    print(len(mutant_df))

    # 3. Calculate attributions for each original smiles
    ids = list(pd.unique(mutant_df['molid']))
    assert len(ids) == len(smiles_list)

    sml_attr = []

    for i in tqdm(range(len(ids))):
        molid = ids[i]
        print(f"Processing {molid}")
        mutant_df_mol = mutant_df[mutant_df['molid'] == molid]
        pred = original_preds[i]
        orig_selfie = pd.unique(mutant_df_mol['selfies'])[0]
        mutant_df_mol['difference'] = pred - mutant_df_mol['predictions']
        mean_of_position = mutant_df_mol.groupby(by=['molid', 'mutation_pos'])["difference"].mean().tolist()
        print('Got selfies attributions.')

        # translate to the SMILES world
        selfie, scores, smiles, smiles_tokens, smiles_scores_list, atom_weights = score_from_selfies_to_smiles(
            orig_selfie, mean_of_position)
        cansmi = Chem.MolToSmiles(MolFromSmiles(smiles), canonical=True)
        order = get_correct_order(smiles)
        ordered_smiles_attributions = np.array(atom_weights)[order]
        print('Got smiles attributions')

        attributes_dict = {"model_prediction": pred, "selfie": selfie, "canonical_smiles": cansmi,
                           "selfies_scores": mean_of_position, "attributions": ordered_smiles_attributions,
                           "original_smiles": smiles_list[i], 'compound_id': molid}

        sml_attr.append(attributes_dict)

    return pd.DataFrame.from_records(sml_attr)
