# Run log: DPD Part B

## Windows: VS Code PowerShell

Choose the actual GCC installed on your machine. A typical MSYS2 UCRT64 path is:

```powershell
$gcc = "C:\msys64\ucrt64\bin\gcc.exe"
& $gcc --version
py --version
py -m pip install numpy matplotlib
```

If your existing course installation uses `mingw64`, substitute its GCC path. If your Python command is `python`, substitute it for `py` in every command. The helper passes compiler arguments directly and handles paths containing spaces.

Run the validation checks first:

```powershell
py scripts/run_b.py checks --cc $gcc --name b_checks
py scripts/run_b.py velocity-check --cc $gcc --name b_velocity_check
```

Produce the equilibrated conservative-on fluid and conservative-off ideal reference. Both use 20,000 steps at dt=0.04, discarding the first 5,000 steps for statistical sampling:

```powershell
py scripts/run_b.py fluid --cc $gcc
py scripts/run_b.py ideal --cc $gcc
```

Run B1 from the same `data/b_fluid.restart` in all four cases. Each reaches physical time 200, with D/R off and C on:

```powershell
py scripts/run_b.py energy --cc $gcc --dt 0.04
py scripts/run_b.py energy --cc $gcc --dt 0.02
py scripts/run_b.py energy --cc $gcc --dt 0.01
py scripts/run_b.py energy --cc $gcc --dt 0.005
```

Create every B plot and result table:

```powershell
py scripts/post_b1.py
py scripts/post_b2.py
py scripts/post_b3.py
```

Alternatively the complete simulation sequence above is:

```powershell
py scripts/run_b.py all --cc $gcc
py scripts/post_b1.py
py scripts/post_b2.py
py scripts/post_b3.py
```

The `all` command runs simulations and checks; it does not automatically run the three post-processing scripts.
