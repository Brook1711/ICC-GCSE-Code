import numpy as np

# calculate NMSE of vec_a and vec_b
def NMSE(vec_a, vec_b):
    # normalize vec_a and vec_b
    # vec_a = vec_a / np.linalg.norm(vec_a)
    # vec_b = vec_b / np.linalg.norm(vec_b)
    
    return np.sum(
        np.power(
            np.abs(vec_a - vec_b), 2
        )
    ) / np.sum(
        np.power(
            np.abs(vec_b), 2
        )
    )
    
def NMSE_v2(vec_a, vec_b):
    # normalize vec_a and vec_b
    vec_a = vec_a / np.linalg.norm(vec_a)
    vec_b = vec_b / np.linalg.norm(vec_b)
    
    return np.sum(
        np.power(
            np.abs(vec_a - vec_b), 2
        )
    ) / np.sum(
        np.power(
            np.abs(vec_b), 2
        )
    )