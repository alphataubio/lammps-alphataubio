
# ----------
# LAMMPS-DFT 
# ----------


## DEPENDENCIES

- libxc
- pseudopod
- kokkoskernels




## NOTES

### `diagonalization='david'` in Quantum Espresso (QE)

In Quantum Espresso (QE), the `diagonalization='david'` parameter specifies the use of the **Davidson algorithm** for diagonalizing the Hamiltonian matrix in electronic structure calculations.

The Davidson algorithm is a popular iterative method for eigenvalue problems, particularly when you need to find a small number of eigenvalues and eigenvectors of a large sparse matrix, which is often the case in electronic structure calculations.

#### Why use `diagonalization='david'`?

1. **Efficiency**: The Davidson algorithm is particularly efficient when you're only interested in the lowest eigenvalues (like in many electronic structure calculations). This is because it allows the computation of the subspace of interest without needing to solve for the entire spectrum of eigenvalues.
   
2. **Memory**: It is more memory-efficient for large systems, as it avoids constructing the full eigenvalue spectrum.

3. **Robustness**: It is known to be robust and can handle challenging cases where other diagonalization methods might struggle.

#### Alternatives

- **`diagonalization='cg'`** (Conjugate Gradient): This method is another iterative technique used in QE for diagonalizing the Hamiltonian, which is sometimes preferred for its own efficiency in different types of problems.

### Conclusion

In short, using `diagonalization='david'` invokes the Davidson method in Quantum Espresso for solving the eigenvalue problem efficiently, especially when dealing with large systems.
```
