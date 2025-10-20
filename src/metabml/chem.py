# TODO: 集成 RDKit/图神经网络特征
def smiles_to_features(smiles: str, dim: int = 128):
    # 占位：返回固定向量
    import numpy as np
    rng = np.random.RandomState(abs(hash(smiles)) % (2**32 - 1))
    return rng.rand(dim).astype("float32")
