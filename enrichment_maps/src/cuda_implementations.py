from numba import cuda
import numpy as np
import math

# CUDA kernel
@cuda.jit
def kernel_func(__flatten, __access, C):
    
    # The 2D threadIdx.x and threadIdx.y are a proxy of the position of the matrix
    row, col = cuda.grid(2)
    if row < C.shape[0] and col < C.shape[1]:
        s_1_idx = __access[row]
        e_1_idx = __access[row+1]
        
        s_2_idx = __access[col]
        e_2_idx = __access[col+1]
        
        # Count the unique values in the sorted arrays s_1 and s_2
        i = s_1_idx
        j = s_2_idx
        shared = 0
        while i < e_1_idx and j < e_2_idx:
            a = __flatten[i]
            b = __flatten[j]
            if a == b:
                shared += 1
                i += 1
                j += 1
            elif a < b:
                i += 1
            else:
                j += 1
        
        # Update the matrix
        C[row, col] = (e_1_idx - s_1_idx) + (e_2_idx - s_2_idx) - shared
        
# Initialize the data arrays
def calculate_overlap_matrix(flattened_bins, accession_idxs, bin_count):
    # Copy the arrays to the device
    flattened_bins_global_mem = cuda.to_device(flattened_bins)
    accession_idxs_global_mem = cuda.to_device(accession_idxs)

    # Allocate memory on the device for the result
    result_global_mem = cuda.device_array((bin_count , bin_count), dtype=np.uint16)

    # Configure the blocks
    threadsperblock = (32, 32)
    blockspergrid_x = int(math.ceil(bin_count / threadsperblock[0]))
    blockspergrid_y = int(math.ceil(bin_count / threadsperblock[0]))
    blockspergrid = (blockspergrid_x, blockspergrid_y)

    # Start the kernel 
    kernel_func[blockspergrid, threadsperblock](flattened_bins_global_mem, accession_idxs_global_mem, result_global_mem)

    # Copy the result back to the host
    return result_global_mem.copy_to_host()