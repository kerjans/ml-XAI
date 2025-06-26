from pathlib import Path

DATA_DIR = Path(__file__).parent.parent/"data"
SOL = Path(DATA_DIR) / 'Lipophilicity.csv'
MODEL = Path(DATA_DIR) / 'lipo_model.pkl'
