import numpy as np
from tqdm import tqdm
import networkx as nx
from tool_box.utils import NMSE, NMSE_v2
# fix the random seed
np.random.seed(0)

# factor_cluster = 2.1 1000, 5db
factor_cluster = 2.1
factor_non_active = 0.005
# factor_non_active_scaling = 0.035 1000, 5db
factor_non_active_scaling = 0.035
factor_active = 1.

class EMRF:
    def __init__(self, Graph : nx.Graph, y, Phi, vec_H_a, xi, l_x_abs_max, l_y_abs_max, Geo_Approx_para, sparseness = 0.1) -> None:
        self.factor_non_active = factor_non_active
        
        # input Graph is a networkx graph
        self.Graph = Graph
        # Graph node set index except the 'a' and 'b' node
        self.xi = xi
        
        # init meta Graph data
        self.l_x_abs_max = l_x_abs_max
        self.l_y_abs_max = l_y_abs_max
        self.Geo_Approx_para = Geo_Approx_para
        
        # OMP init
        self.y = y
        self.Phi = Phi
        self.hat_h_sl = np.mat(
            np.zeros(
                (len(xi),1)
            ), dtype= np.complex128
        )
        
        # init residual
        self.res = y
        
        # init temperary target
        self.x = np.reshape(
            np.zeros(Phi.shape[1], dtype=np.complex128),
            (-1, 1)
        )
        
        # init pseudo matrix
        self.Phi_pinv = np.linalg.pinv(self.Phi)
        self.Phi_H = self.Phi.H
        
        # init Geo_Approx_para
        self.vec_H_a = vec_H_a
        vec_H_a_power = np.sum(
            np.power(
                np.abs(vec_H_a), 2
            )
        )
        self.Geo_Approx_para['variance_beta'] = self.factor_non_active
        self.K_true, self.tau_true, self.alpha_true = self.find_real_sparseness(vec_H_a)
        self.Geo_Approx_para['variance_alpha'] = vec_H_a_power / self.K_true * factor_active
        
        # result container
        self.NMSE_list = []
        self.NMSE_list_v2 = []
        
        pass
    
    def alpha_beta_swap(self,):
        # 0.1 update residue
        self.update_residue()
        
        # update self.factor_non_active according to the residue
        res_power = np.sum(
            np.power(
                np.abs(self.res), 2
            )
        )   / len(self.res)
        self.factor_non_active = res_power * factor_non_active_scaling
        
        # 0.2 update temperary target
        self.update_temperary_target()
        
        # 1.1 update Geo_Approx_para
        # self.update_Geo_Approx_para()
        
        # 1.2 update edge weights
        self.update_Graph()
        
        # 2.1 find the minimum_cut
        self.min_cut_set = self.get_mini_cut_set()
        
        # 2.2 update the support
        self.update_support()
        
        # 0.3 update estimated h_sl
        self.update_h_sl()
        
        self.NMSE_list_v2.append(
            NMSE_v2(
                np.reshape(self.vec_H_a, (-1,1)), 
                np.reshape(self.hat_h_sl, (-1,1))
            )
        )
        self.NMSE_list.append(
            NMSE(
                np.reshape(self.vec_H_a, (-1,1)),
                np.reshape(self.hat_h_sl, (-1,1))
            )
        )
        
        # self.factor_non_active = self.factor_non_active * 0.9
        
        return self.hat_h_sl
    
    # update Geo_Approx_para
    def update_Geo_Approx_para(self,):
        K, tau, alpha = self.find_real_sparseness(self.x)
        self.Geo_Approx_para['K'] = K
        self.Geo_Approx_para['tau'] = tau
        self.Geo_Approx_para['alpha'] = alpha
        return
    
    # update Graph
    def update_Graph(self,):
        # update Graph edge weight between 'a', 'b' nodes and l-th node
        for l_idx, l_pair in enumerate(self.xi):
            node_idx = self.from_l_to_idx(l_pair)
            node = self.Graph.nodes[node_idx]
            
            ##  update the edges between the 'a' and l-th node
            self.Graph['a'][node_idx]['weight'] = \
                self.func_D_v2(
                    h_sl = self.x[l_idx, 0],
                    s_l = 1
                ) + \
                self.func_V(
                    s_l = 1,
                    s_l_prime = node['s']
                )
            
            ##  update the edges between the 'b' and l-th node
            self.Graph['b'][node_idx]['weight'] = \
                self.func_D_v2(
                    h_sl = self.x[l_idx, 0],
                    s_l = -1
                ) + \
                self.func_V(
                    s_l = -1,
                    s_l_prime = node['s']
                )
        # update Graph edge weight between l-th node and l'-th node
        for l_idx, l_pair in enumerate(self.xi):
            node_idx = self.from_l_to_idx(l_pair)
            node = self.Graph.nodes[node_idx]

            for nei_node in self.neighbor_l(l_pair):
                nei_node_idx = self.from_l_to_idx(nei_node)
                nei_node = self.Graph.nodes[nei_node_idx]
                self.Graph[node_idx][nei_node_idx]['weight'] = \
                    self.func_V(
                        s_l = node['s'],
                        s_l_prime = nei_node['s']
                    )
        
        return
    
    # update res by OMP
    def update_residue(self,):
        self.res = self.y - np.reshape(self.Phi @ self.hat_h_sl, (-1, 1))
        return
    
    # update x
    def update_temperary_target(self,):
        # self.x = np.reshape(
        #     self.Phi_pinv @ self.res,
        #     (-1,1)
        # ) + self.hat_h_sl
        
        self.x = np.reshape(
            self.Phi_H @ self.res,
            (-1,1)
        ) + self.hat_h_sl
        return
    
    # update support
    def update_support(self,):
        # update support variable
        for edge in self.min_cut_set:
            if 'a' not in edge and 'b' not in edge:
                continue
            
            if 'a' in edge:
                for node_idx in edge:
                    if node_idx != 'a':
                        self.Graph.nodes[node_idx]['s'] = 1
                continue
            
            if 'b' in edge:
                for node_idx in edge:
                    if node_idx != 'b':
                        self.Graph.nodes[node_idx]['s'] = -1

        return
    
    # update h_sl
    def update_h_sl(self,):
        support = self.get_support()
        activated_idx = np.where(support == 1)[0] 
        
        self.hat_h_sl[activated_idx] = np.linalg.pinv(
            self.Phi[:, activated_idx]
        ) @ self.y
        
        # truncate the h_sl by choosing the largest K elements
        hat_h_sl_abs = np.abs(self.hat_h_sl)
        self.hat_h_sl[np.where(hat_h_sl_abs < self.tau_true)[0]] = 0
        
        return self.hat_h_sl
    
    # get
    def get_support(self,):
        # init support array
        support = np.zeros(
            (len(self.xi),1)
        )
        
        # get support
        for l_idx, l_pair in enumerate(self.xi):
            node_idx = self.from_l_to_idx(l_pair)
            node = self.Graph.nodes[node_idx]
            support[l_idx] = node['s']
        
        return support
    
    # Energy Function

    def func_D_v2(self, h_sl, s_l):
        variance_alpha = self.Geo_Approx_para['variance_alpha']
        sigma_alpha = np.sqrt(variance_alpha)
        # variance_beta = self.Geo_Approx_para['variance_beta']
        variance_beta = self.factor_non_active
        sigma_beta = np.sqrt(variance_beta)
        h_sl_abs = np.abs(h_sl)
        
        if s_l == 1:
            # p(h_sl | s_l = 1) is a complex Gaussian distribution with 0 mean and variance_alpha variance
            pdf = (1/(np.pi * variance_alpha)) * np.exp(- (h_sl_abs ** 2) / variance_alpha)
            return - np.log10(pdf)
        elif s_l == -1:
            # p(h_sl | s_l = -1) is a complex Gaussian distribution with 0 mean and variance_beta variance
            pdf = (1/(np.pi * variance_beta)) * np.exp(- (h_sl_abs ** 2) / variance_beta)
            return - np.log10(pdf)
        else:
            raise ValueError('s_l must be 1 or -1')

    
    def func_D(self, h_sl, s_l):
        Geo_Approx_para = self.Geo_Approx_para
        tau = Geo_Approx_para['tau']
        if s_l == 1:
            if h_sl <= tau:
                epsilon_4 = Geo_Approx_para['epsilon_4']
                return - np.log10(epsilon_4)
            else:
                epsilon_3 = Geo_Approx_para['epsilon_3']
                return - np.log10(epsilon_3)
        elif s_l == -1:
            if h_sl <= tau:
                epsilon_1 = Geo_Approx_para['epsilon_1']
                return - np.log10(epsilon_1)
            else:
                epsilon_2 = Geo_Approx_para['epsilon_2']
                return - np.log10(epsilon_2)
        else:
            raise ValueError('s_l must be 1 or -1')
        
    def func_V(self, s_l, s_l_prime, eta = 0.15 * factor_cluster):
        return eta * (2.1 - s_l * s_l_prime)

    # Util Functions
    def from_idx_to_l(self, idx):
        return (idx[0] - self.l_x_abs_max, idx[1] - self.l_y_abs_max)

    def from_l_to_idx(self, l):
        return (l[0] + self.l_x_abs_max, l[1] + self.l_y_abs_max)
    
    def neighbor_l(self, l):
        Graph = self.Graph
        (l_x_abs_max, l_y_abs_max) = (self.l_x_abs_max, self.l_y_abs_max)
        (l_x, l_y) = l
        node_idx = (l_x + l_x_abs_max, l_y + l_y_abs_max)
        res = []
        for nei_node_idx in Graph.neighbors(node_idx):
            if nei_node_idx == 'a' or nei_node_idx == 'b':
                continue
            else:
                res.append(Graph.nodes[nei_node_idx]['l'])
        return res
    
    # find minimum cut set
    def get_mini_cut_set(self,):
        cut_value, partition = nx.minimum_cut(self.Graph, "a", "b", capacity = 'weight')
        reachable, non_reachable = partition

        cutset = set()
        for u, nbrs in ((n, self.Graph[n]) for n in reachable):
            cutset.update((u, v) for v in nbrs if v in non_reachable)
        return cutset
    
    def find_real_sparseness(self, vec_H_a, threshold=0.9):
        vec_H_a_abs = np.abs(vec_H_a)
        M = vec_H_a_abs.shape[0]
        # sort the vector in descending order
        vec_H_a_sorted = np.sort(np.reshape(np.array(vec_H_a_abs), (M,)), )
        vec_H_a_sorted = np.flip(vec_H_a_sorted, axis=0)
        overall_power = np.sum(vec_H_a_sorted)
        
        K = 0
        for i, val in enumerate(vec_H_a_sorted):
            if np.sum(vec_H_a_sorted[0:i]) >= threshold * overall_power:
                K = i
                break
            
        if 2 * K > M:
            tau = vec_H_a_sorted[-10]
        else:
            tau = vec_H_a_sorted[K]
        alpha = vec_H_a_sorted[0]
        return K, tau, alpha
    
    def get_NMSE(self,):
        return NMSE_v2(self.Phi @ self.vec_H_a, self.Phi @ self.hat_h_sl)