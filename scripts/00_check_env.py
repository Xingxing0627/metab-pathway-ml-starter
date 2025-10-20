import sys, json, importlib, platform

pkgs = ["numpy","pandas","sklearn","torch","rdkit","yaml","tqdm","matplotlib","networkx","faiss"]
info = {"python": sys.version, "platform": platform.platform(), "packages": {}}
for p in pkgs:
    try:
        m = importlib.import_module(p if p!="yaml" else "yaml")
        if p=="torch":
            cuda = getattr(m, "cuda", None)
            info["packages"][p] = {"version": getattr(m, "__version__", "unknown"),
                                   "cuda_available": bool(cuda and cuda.is_available())}
        else:
            info["packages"][p] = {"version": getattr(m, "__version__", "unknown")}
    except Exception as e:
        info["packages"][p] = {"error": str(e)}
print(json.dumps(info, indent=2))
