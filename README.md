# IsoSolve notebook
IsoSolve is a computation framework to improve isotopic coverage and consolidate isotopic measurements by MS and/or NMR. IsoSolve is available from [this repository](https://github.com/MetaSys-LISBP/IsoSolve).

IsoSolve can be used as a Python module that you can import directly, for instance in [Jupyter notebooks](https://test-jupyter.readthedocs.io/en/latest/install.html) or in your own software.
We showcase IsoSolve usage with step-by-step examples in a Jupyter notebook distributed via this repository. We also distribute an HTML file showing the [notebook’s output after execution](https://htmlpreview.github.io/?https://github.com/MetaSys-LISBP/IsoSolve_notebook/blob/main/html/IsoSolve_notebook.html).

# Installation and usage

If not yet done, you can install Jupyter with:

```bash
pip install --user jupyter
```

Some dependencies are specific to our notebook, not to IsoSolve itself. If these dependencies are not already available on your system, they can be installed with:

```bash
pip install --user isosolve seaborn matplotlib
```

Download and unpack the notebook's [tarball](https://github.com/MetaSys-LISBP/IsoSolve_notebook/archive/main.tar.gz) and go in a shell to the notebook's directory.

After that, you are ready to examine and execute the notebook by launching:

```bash
jupyter notebook IsoSolve.ipynb
```

The notebook will open in your web browser where in each cell you can read/modify/execute a proposed code as well as read accompanying comments.
