def pick_device(preferred: str | None = None) -> str:
    """
    Detecta el dispositivo disponible para PyTorch.
    Devuelve 'cuda', 'mps' (Apple) o 'cpu'.
    Si se pasa preferred='cpu'/'cuda'/'mps', lo fuerza.
    """
    if preferred in {"cpu", "cuda", "mps"}:
        return preferred

    try:
        import torch
    except Exception:
        return "cpu"

    try:
        if torch.cuda.is_available():
            return "cuda"
    except Exception:
        pass

    try:
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
    except Exception:
        pass

    return "cpu"
