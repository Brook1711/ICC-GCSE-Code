# Orthogonal Matching Pursuit algorithm
import numpy as np
# tqdm
from tqdm import tqdm
from tool_box.utils import NMSE, NMSE_v2
# fix the seed  
np.random.seed(0)

class OMP:
    def __init__(self, A, y, vec_H_a, sparseness = 0.5, norm_err_threshold = 1e-1):
        self.A = A  # A is the transform matrix
        self.y = y  # y is the observed signal
        self.k = int(sparseness * A.shape[1])  # k is the sparsity
        self.x = np.zeros(A.shape[1], dtype=np.complex128)
        self.res = y # res is the residual
        self.idx = [] # idx is the index of the selected atoms
        self.err = [] # err is the error
        y_norm = np.linalg.norm(y)
        self.vec_H_a = vec_H_a
        self.err_threshold = norm_err_threshold * y_norm
        self.iter = 0 # iter is the number of iterations
        self.pinv_A = np.linalg.pinv(self.A)
        
        self.NMSE_list = []
        self.NMSE_list_v2 = []

    def run(self):
        for i in tqdm(range(self.k)):
            self.idx.append(
                np.argmax(np.abs(
                    self.pinv_A @ self.res)
                )
            )
            self.x[self.idx] = np.reshape(np.linalg.pinv(self.A[:, self.idx]) @ self.y, (len(self.idx),))
            self.res = self.y - np.reshape(self.A @ self.x, (-1, 1))
            self.err.append(np.linalg.norm(self.res))
            # if self.err[-1] < self.err_threshold:
            #     break
            self.NMSE_list_v2.append(
                NMSE_v2(
                    np.reshape(self.vec_H_a, (-1,1)), 
                    np.reshape(self.x, (-1,1))
                )
            )
            self.NMSE_list.append(
                NMSE(
                    np.reshape(self.vec_H_a, (-1,1)),
                    np.reshape(self.x, (-1,1))
                )
            )
            
            
        return self.x, self.idx, self.err
        # while self.iter < self.k:
        #     self.iter += 1
        #     self.idx.append(
        #         np.argmax(np.abs(
        #             np.linalg.pinv(self.A) @ self.res)
        #         )
        #     )
        #     self.x[self.idx] = np.linalg.pinv(self.A[:, self.idx]) @ self.y
        #     self.res = self.y - np.reshape(self.A @ self.x, (-1, 1))
        #     self.err.append(np.linalg.norm(self.res))
        #     if self.err[-1] < self.err_threshold:
        #         break
        # return self.x, self.idx, self.err

    def get_x(self):
        return self.x

    def get_idx(self):
        return self.idx

    def get_err(self):
        return self.err