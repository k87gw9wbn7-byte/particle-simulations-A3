# Part C run log: Windows VS Code and Jupyter


## 1. Checks

```powershell
py scripts/run_c.py checks --cc $gcc
py scripts/run_c.py tensor-check --cc $gcc
```

Wait for each `Finished` message. Detailed C output is saved in `data/checks_run.log` and `data/tensor-check_run.log`. The helper recompiles automatically for every case. Compiler warnings are errors intentionally. The supplied `%zu` format fixes the known Windows `%lu` issue.

## 2. Production run settings

The helper supplies settings to C at compilation; you do not edit them between runs. All cases use dt=0.04, steps=120000, equil=20000, sample=50 and block=5000. Both C and D/R stay enabled. C1 uses 3000 beads in a cubic 10 box; the other cases use 6000 in 20×10×10 with a separated start. Profile width is 0.25. aAA=aBB=25 throughout. No input restart is loaded.

Run C1:

```powershell
py scripts/run_c.py c1 --cc $gcc
```

Run C2:

```powershell
py scripts/run_c.py c2 --cc $gcc --aab 37
```

Run the four C3/C4 cases:

```powershell
py scripts/run_c.py c3 --cc $gcc --aab 35
py scripts/run_c.py c3 --cc $gcc --aab 39
py scripts/run_c.py c3 --cc $gcc --aab 43
py scripts/run_c.py c3 --cc $gcc --aab 47
```

Run the C4 one-phase control:

```powershell
py scripts/run_c.py control --cc $gcc
```

Each command prints compile/run/finished messages. C progress is redirected to `data/<prefix>_run.log`; a quiet terminal while it runs is expected. For progress in a second PowerShell terminal, for example:

```powershell
Get-Content data/c2_run.log -Tail 8
```

To run all checks and all seven production cases sequentially:

```powershell
py scripts/run_c.py all --cc $gcc
```

C3 and C4 use the same four demixed runs. Do not run a duplicate four-case C4 series. The control is the only extra C4 simulation.

## 3. Notebook setup and execution

Open A3.ipynb and select your Python kernel. If needed, install into that kernel:

```python
%pip install numpy matplotlib
```

Check the working directory:

```python
from pathlib import Path
print(Path.cwd())
assert Path('data').is_dir() and Path('code').is_dir()
```

If necessary:

```python
import os
os.chdir(r'C:\Users\20224045\OneDrive - TU Eindhoven\Documents\Master\PBS\A3\particle-simulations-A3')
```

Copy the complete contents of each notebook-ready Python file into its corresponding C1/C2/C3/C4 code cells. The files are self-contained, use `Path.cwd()`, have no `__file__`, and call `plt.show()`, so plots appear inline. Do not put Python code inside a Markdown cell or paste surrounding Markdown code fences into a code cell.

Alternatively, keep those files and run notebook cells:

```python
%run notebook_cells/C1_RDF.py
```

```python
%run notebook_cells/C2_demixing.py
```

```python
%run notebook_cells/C3_chi.py
```

```python
%run notebook_cells/C4_surface_tension.py
```

Run C3 before C4: C4 reads the chi estimates and per-run discard times saved by C3. The four `.py` files are designed for notebooks, not simulation launchers. Running them does not change C settings or rerun a simulation.

