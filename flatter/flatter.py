def flatter(M):
    import os
    import re
    from subprocess import check_output
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(os.cpu_count())
    ret = subprocess.check_output(["flatter"], input=z.encode(), env=env)
    return matrix(M.nrows(), M.ncols(), map(int, findall(rb"-?\d+", ret)))
