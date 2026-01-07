
# Github Supplementary Material (Anonymous)

This repository provides supplementary materials for an **anonymous submission under review**, including the experimental code and reference circuit implementations used in the paper.



## Requirements 
(Tested on a multi-core system with 32 CPU cores and 128 GB of RAM, running Ubuntu 24.04 LTS.)

* C++17-compatible compiler (e.g., `g++ >= 9`)
* Gurobi Optimizer (C++ API, tested with version 12.0): https://www.gurobi.com
* Standard Unix utilities (`bash`, `mkdir`, `rm`, `tee`)



## Code (`./code`)

### Usage

```bash
sh run.sh <OPS_FILE> <DEPTH_LIMIT> <TIE_BREAKING> <TARGET_STRATEGY>
```

### Arguments

* **OPS_FILE**: S-box instance

  ```
  aes10_ops, aes10v2_ops, aes10v3_ops
  aes12_ops, aes12v2_ops, aes12v3_ops
  snow_ops, saturnin_ops
  ```

In `./code/data/*_ops.txt`, each line has the format `<id src1 src2 type flag>`,
where `type in {xor, nxor, and, or}`, and `flag = -1` for non-output nodes; otherwise, `flag` denotes the output index.

* **DEPTH_LIMIT**: Maximum circuit depth

* **TIE_BREAKING** (Sec. 5.2):

  * `1` MN (Maximum Norm)
  * `2` ASM (Average Support Minimization)
  * `3` SOM (Support Overlap Maximization)

* **TARGET_STRATEGY** (Sec. 5.2):

  * `1` AT (All targets)
  * `2` ST (Selective targets)

### Example

```bash
sh run.sh aes10_ops 18 1 2
```


### Directory Structure

```
code/
├── src/                        # C++ implementation
│   ├── main.cpp                # Program entry
│   ├── super_bp.inl            # Core BP heuristics 
│   ├── my_refine.inl           # Local refinement algorithm
│   ├── parallel_algorithms.hpp # Parameter control
│   ├── thread_pool.hpp         # Thread pool for multi-core execution
│   └── *.hpp / *.inl           # Utility headers and helpers
│
├── s-box/                      # S-box functional descriptions
│   ├── aes.txt
│   ├── snow3g.txt
│   └── saturn.txt
│
├── data/                       # Converted S-box circuits (OPS format) circuits
│   └── *_ops.txt 
│
└── run.sh                      # Experiment driver (compile and run)
```




## Code for Runtime Comparison (`./others`)

The runtime comparison in the rebuttal is based on three implementations code:

1. **Original baseline ([JBKK25], Python)**  
The original implementation is available at  https://github.com/lemontrr/Extended-BP-Framework  and was used without modification (mirrored in `others/baseline_jbkk25/` for convenience).

2. **Baseline reimplementation (C++)**  
A C++ reimplementation of [JBKK25], independently implemented by us (provided in `others/baseline_jbkk25_cpp_reimp/`), used to evaluate the effect of language-level optimization. **Usage:** `sh run.sh <OPS_FILE> <DEPTH_LIMIT>`.

3. **Improved method (this work)**  
The improved C++ implementation proposed in this paper (provided in `code/`), used to obtain the reported speed-up results.



## Circuit Implementations (`./circuit_implementations`)

Python implementations of S-box circuits achieving the **minimum XOR counts** reported in the paper (AES, SNOW3G, Saturnin).


### 1. Included S-Box Implementations: 
The S-Box circuit implementations used in our paper:

AES S-box
- **AES10** — Circuit from [BP10]
- **AES12** — Circuit from [BP12]
- **AES10v2**, **AES10v3** — Circuits from [JBKK25]
- **AES12v2**, **AES12v3** — Circuits from [JBKK25]

Saturnin Super S-box
- **Saturnin** — Circuit from [CDL+20]

SNOW3G S-box
- **SNOW** — Circuits from [JBKK25], [BHNS10]



### 2. Circuit Implementation

All Python files follow a unified structure. The file names encode the S-box type, depth bound (H), and XOR count. For example:

```
AES10_H18_XORS91.py means: AES10 circuit, depth bound H = 18, XOR gate count = 91.
```
Each file contains two main components:

#### 2.1 S-Box Circuit Implementation  
Functions such as `AES_Sbox()`, `Saturnin_Super_Sbox()`, or `SNOW_SQ()` perform the Boolean circuit computation.

##### Input format
Each S-box function processes a single **8-bit** input:

```python
x = [(X >> i) & 1 for i in range(8)].
```
The input byte is decomposed into 8 individual bits for XOR/AND Boolean operations.

##### Intermediate variables

A sequence of temporary nodes (e.g., t[i], r[i]) represents the Boolean gates.
Each line corresponds to one XOR or AND operation, exactly following the structure given in our paper.

##### Output format

The 8 output bits are combined back into an 8-bit integer:

```python
Y = (y[7]<<7) | (y[6]<<6) | ... | (y[0]<<0).
```
producing the final S-box output value (0–255).

#### 2.2 Correctness Verification

Each file includes a verification script (adapted from [JBKK25]) that:

+ Computes all S-box outputs

+ Compares them against the official reference implementation

+ Outputs "Right!!" or "Wrong!!"

+ Reports detailed mismatches if discrepancies are found

This ensures the Boolean circuit implementations are functionally correct.
