# Baseline reimplementation (C++)

This is a C++ reimplementation of [JBKK25]

## Requirements
Dependencies

+ C++17-compatible compiler (e.g., g++ >= 9)
+ Gurobi Optimizer (C++ API, tested with version 12.0): https://www.gurobi.com
+ Standard Unix utilities (bash, mkdir, rm, tee)


## Usage

```bash
sh run.sh <filename> <depth>
```
## Parameters

+ filename: Input S-box file name, from `code_target_imps` directory. Available options: 
```
AES10, AES10v2, AES10v3, ...
```
+ depth: Circuit depth constraints

## Example
```bash
sh run.sh AES10 23
```

## Output Files

+ Log Files
Log files are saved in the `log/` directory

+ Result Files
Generated results are saved in the `Results/` directory.
