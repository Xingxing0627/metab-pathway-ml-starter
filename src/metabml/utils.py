def set_seed(seed: int = 42):
    import random, numpy as np
    random.seed(seed); np.random.seed(seed)
    try:
        import torch; torch.manual_seed(seed)
        if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)
    except Exception:
        pass
