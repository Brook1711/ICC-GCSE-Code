# %%
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from tqdm import tqdm
from tool_box.utils import NMSE
# fix random seed
np.random.seed(0)

plt.rcParams['figure.dpi'] = 300
N_RF = 1000
SNR = 50 # dB

# %% [markdown]
# #### Read Generated Channel Data

# %%
# read data from \generate_channel_data\generate_channel_data.mat
import scipy.io as sio
generate_channel_data = sio.loadmat('./data/generate_channel_data.mat')
cluster_para = generate_channel_data['meta_data']['cluster_para']
# xi = generate_channel_data['meta_data']['xi']
f_c = generate_channel_data['meta_data']['f_c']
N_x = int(generate_channel_data['meta_data']['N_x'])
N_y = int(generate_channel_data['meta_data']['N_y'])
lambda_c = generate_channel_data['meta_data']['lambda_c']
delta = generate_channel_data['meta_data']['delta']
N = int(generate_channel_data['meta_data']['N'])
L_x = generate_channel_data['meta_data']['L_x']
L_y = generate_channel_data['meta_data']['L_y']

H_channel = np.mat(generate_channel_data['channel']['H_channel'][0][0])
vec_H_a = np.mat(generate_channel_data['channel']['vec_H_a'][0][0])
vec_H_AD = np.mat(generate_channel_data['channel']['vec_H_AD'][0][0])
variance = np.mat(generate_channel_data['channel']['variance'][0][0])
# Psi = np.mat(generate_channel_data['channel']['Psi'][0][0])
# Psi_AD = np.mat(generate_channel_data['channel']['Psi_AD'][0][0])

K = 0
sparseness = 0.1


sparse_power_threshold = 0.9
eta = 0.15
Geo_Approx_para = {
    'epsilon_1' : 1e4,
    'epsilon_2' : 1e-3,
    'epsilon_3' : 1e3,
    'epsilon_4' : 1e-2,
    'K' : 0,   # the number of nonZero elements in the sparse vector
    'tau' : 0.1,    # the 2K largest entries of the sparse vector
    'alpha' : 5     # the biggest entry of the sparse vector
}

# %%
xi = []
l_x_abs_max = int(L_x / lambda_c)
l_y_abs_max = int(L_y / lambda_c)
print("l_x_abs_max = ", l_x_abs_max)
print("l_y_abs_max = ", l_y_abs_max)

# for start with -l_x_abs_max to l_x_abs_max
for l_x in range(-l_x_abs_max, l_x_abs_max + 1):
    # for start with -l_y_abs_max to l_y_abs_max
    for l_y in range(-l_y_abs_max, l_y_abs_max + 1):
        k_x = 2 * np.pi * l_x / L_x
        k_y = 2 * np.pi * l_y / L_y
        if k_x ** 2 + k_y ** 2 < (2 * np.pi / lambda_c) ** 2:
            xi.append((l_x, l_y))


xi_AD = [] # angular domain index pair

n_x_abs_max = int((N_x - 1) / 2)
n_y_abs_max = int((N_y - 1) / 2)
print("n_x_abs_max = ", n_x_abs_max)
print("n_y_abs_max = ", n_y_abs_max)

# for start with -n_x_abs_max to n_x_abs_max
for n_x in range(-n_x_abs_max, n_x_abs_max + 1): 
    # for start with -n_y_abs_max to n_y_abs_max
    for n_y in range(-n_y_abs_max, n_y_abs_max + 1):
        xi_AD.append((n_x, n_y))

M = len(xi_AD)
print("M = ", M)

# %% [markdown]
# #### Regenerate `Psi_AD`

# %%
def Psi_AD_col(n_prime_pair):
    n_prime_x = n_prime_pair[0]
    n_prime_y = n_prime_pair[1]
    amp = 1 / np.sqrt(N)
    exp_func = np.exp(
        1j * (
            2*np.pi * n_prime_x / N_x * (np.arange(N_x) - (N_x - 1) / 2) + \
                2*np.pi * n_prime_y / N_y * (np.arange(N_y) - (N_y - 1) / 2)[:, np.newaxis]
        )
    )
    vec = np.reshape(exp_func, (N_x*N_y, ))
    return amp * vec

xi_AD_norm = len(xi_AD)
Psi_AD = np.zeros((N, xi_AD_norm), dtype=np.complex_)
for n_prime, n_prime_pair in enumerate(xi_AD):
    Psi_AD[:, n_prime] = Psi_AD_col(n_prime_pair)
    # print(n_prime_pair)
Psi_AD = np.mat(Psi_AD)

# %% [markdown]
# #### Visualize Channel Function

# %%
# visualize the channel
def visualize_wavenumber_domain_channel(vec_H_a):
    mat_H_a = np.mat(np.zeros((2*l_x_abs_max+1, 2*l_y_abs_max+1), dtype=complex))
    # all set to np.nan in mat_H_a
    mat_H_a[:,:] = np.nan
    for l, l_pair in enumerate(xi):
        l_x = l_pair[0]
        l_y = l_pair[1]
        mat_H_a[l_x + l_x_abs_max, l_y + l_y_abs_max] = vec_H_a[l]

    # plt.figure(figsize=(10, 10))
    plt.imshow(np.abs(mat_H_a), extent=[-l_x_abs_max, l_x_abs_max, -l_y_abs_max, l_y_abs_max], cmap='hot')
    plt.colorbar()
    plt.title('Instantaneous channel in WD')
    
def visualize_wavenumber_variance(variance):
    # plt.figure(figsize=(8, 6))
    plt.imshow(variance, cmap='jet', extent=[-l_x_abs_max, l_x_abs_max, -l_y_abs_max, l_y_abs_max])
    plt.colorbar()
    plt.title('Variance of the channel')

def DFT_power_spectrum(H_channel):
    # DFT power spectrum of the channel
    DFT_power = np.abs(np.fft.fft2(np.reshape(H_channel, (N_x, N_y))))

    # element shift according to the DFT power spectrum
    DFT_power = np.fft.fftshift(DFT_power)
    
    # normalize the DFT power spectrum
    DFT_power = DFT_power / np.max(DFT_power)

    # visualize the DFT power spectrum
    # plt.figure(figsize=(10, 10))
    plt.imshow(DFT_power, cmap='jet', extent=[-np.pi/2, np.pi/2, -np.pi/2, np.pi/2])
    plt.colorbar()
    plt.title('normmalized DFT power spectrum')
    return DFT_power

# visualize the angular domain channel
def visualize_angular_domain_channel(vec_H_AD, trim = False):
    mat_H_AD = np.mat(np.zeros((N_x, N_y), dtype=complex))
    # all set to np.nan in mat_H_a
    # mat_H_a[:,:] = np.nan
    for n_prime, n_prime_pair in enumerate(xi_AD):
        n_prime_x = n_prime_pair[0]
        n_prime_y = n_prime_pair[1]
        mat_H_AD[n_prime_x + n_x_abs_max, n_prime_y + n_y_abs_max] = vec_H_AD[n_prime]

    # plt.figure(figsize=(10, 10))
    plt.imshow(np.abs(mat_H_AD), cmap='jet', extent=[- n_x_abs_max, n_x_abs_max, - n_y_abs_max, n_y_abs_max])
    
    # trim the image only show the center part
    if trim:
        plt.xlim(-l_x_abs_max, l_x_abs_max)
        plt.ylim(-l_y_abs_max, l_y_abs_max)
    plt.colorbar()
    
def find_real_sparseness(vec_H_a, threshold=0.9):
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
    tau = 0
    if 2 * K >= M:
        tau = vec_H_a_sorted[-1]
    else:
        tau = vec_H_a_sorted[2 * K]
    alpha = vec_H_a_sorted[0]
    sparseness = K / M
    return sparseness

# %% [markdown]
# #### Visualize Channel

# %%
# 1 row 3 columns
plt.figure(figsize=(20, 4))
# sub out of 3 in one row
plt.subplot(141)
visualize_wavenumber_domain_channel(vec_H_a)

plt.subplot(142)
visualize_wavenumber_variance(variance)

plt.subplot(143)
visualize_angular_domain_channel(vec_H_AD)

plt.subplot(144)
visualize_angular_domain_channel(vec_H_AD, trim=True)

# print sparseness
print("sparseness = ", find_real_sparseness(vec_H_a))
print("sparseness = ", find_real_sparseness(vec_H_AD))

# %% [markdown]
# #### Observation Model

# %% [markdown]
# $$
# {\bf y} = {\bf C} {\bf \Psi} {\bf h}_a + {\bf n}
# $$

# %%
# Measurement matrxi for compressed sensing, {N}_{RF} \times N
measurement_matrix = np.random.randn(N_RF, N) * 1 / np.sqrt(N)

y = np.dot(
    measurement_matrix,
    np.dot(
        Psi_AD, vec_H_AD
    )
)
Phi_AD = np.dot(measurement_matrix, Psi_AD)

# %% [markdown]
# #### Index-Wavenumber Remapping

# %%
def from_idx_to_l(idx):
    return (idx[0] - l_x_abs_max, idx[1] - l_y_abs_max)

def from_l_to_idx(l, l_x_abs_max = l_x_abs_max, l_y_abs_max = l_y_abs_max):
    return (l[0] + l_x_abs_max, l[1] + l_y_abs_max)

def from_n_prime_to_idx(n_prime):
    return (n_prime[0] + n_x_abs_max, n_prime[1] + n_y_abs_max)

# %% [markdown]
# #### Graph Edge Value Calculation

# %% [markdown]
# $$
# \eta_{l,l^\prime} = 0.5 \arcsin \left( \left( 1-m^8 \right)^{-0.25} \right)
# $$

# %% [markdown]
# $$
# m = \frac{(+1) \times \tilde{K} + (-1) \times (N-\tilde{K})}
# {N}
# $$

# %%
# Useless function
def get_eta(sparseness):
    m = 1 * sparseness + (-1) * (1 - sparseness)
    eta = 0.5 * np.arcsin(
        np.power(
            1 - np.power(m, 8), -0.25
        )
    )
    return eta
# eta = get_eta(sparseness = 0.2)

# %% [markdown]
# #### Graph 'Energy Function' Model

# %% [markdown]
# $$
# \begin{aligned}
# V_{l,l^\prime}\left(s_l, s_{l^\prime}\right) &= - \eta_{l,l^\prime} \cdot \left(s_l s_{l^\prime} - 1 \right) \\
# &= \eta_{l,l^\prime} \cdot \left(1 - s_l s_{l^\prime} \right)
# \end{aligned}
# $$

# %%
def func_V(s_l, s_l_prime):
    return eta * (1 - s_l * s_l_prime)

# %% [markdown]
# $$
# D_l\left( s_l \right) = -\log \left( \hat{p} \left( h^{(j-1)}_{s,l} \right) \right)
# $$

# %%
def func_D(h_sl, s_l, Geo_Approx_para):
    tau = Geo_Approx_para['tau']
    if s_l == 1:
        if h_sl <= tau:
            epsilon_4 = Geo_Approx_para['epsilon_4']
            return np.log10(epsilon_4)
        else:
            epsilon_3 = Geo_Approx_para['epsilon_3']
            return np.log10(epsilon_3)
    elif s_l == -1:
        if h_sl <= tau:
            epsilon_1 = Geo_Approx_para['epsilon_1']
            return np.log10(epsilon_1)
        else:
            epsilon_2 = Geo_Approx_para['epsilon_2']
            return np.log10(epsilon_2)
    else:
        raise ValueError('s_l must be 1 or -1')

# %% [markdown]
# #### Graph init

# %%
# myGraph = nx.grid_2d_graph(int(2*l_x_abs_max+1), int(2*l_y_abs_max+1))

myGraph_AD = nx.grid_2d_graph(int(2*n_x_abs_max+1), int(2*n_y_abs_max+1))

for n_prime_x in range(-n_x_abs_max, n_x_abs_max+1):
    for n_prime_y in range(-n_y_abs_max, n_y_abs_max+1):
        node_idx = from_n_prime_to_idx((n_prime_x, n_prime_y))
        myGraph_AD.nodes[node_idx]['l'] = (n_prime_x, n_prime_y)
        myGraph_AD.nodes[node_idx]['s'] = -1


# add node 'a' and 'b'
myGraph_AD.add_node('a')   # node alpha
myGraph_AD.add_node('b')   # node beta

# add edge between node and 'a', 'b'
for n_prime in xi_AD:
    node_idx = from_n_prime_to_idx(n_prime)
    myGraph_AD.add_edge('a', node_idx, weight = 0)
    myGraph_AD.add_edge('b', node_idx, weight = 0)

# %% [markdown]
# #### Find the non Zero enries in the Sparse Channel

# %%
def find_real_sparseness(vec_H_a, threshold=0.9):
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
    tau = 0
    if 2 * K >= M:
        tau = vec_H_a_sorted[-1]
    else:
        tau = vec_H_a_sorted[2 * K]
    alpha = vec_H_a_sorted[0]
    return K, tau, alpha

K, tau, alpha = find_real_sparseness(vec_H_a, threshold=sparse_power_threshold)
sparseness = K / M
print(
    "K = {}, \nN = {}, \n2K = {} \ntau = {}, \nalpha = {}, \nsparseness = {}".format(K, N, 2*K, tau, alpha, sparseness)
)
Geo_Approx_para['K'] = K
Geo_Approx_para['tau'] = tau
Geo_Approx_para['alpha'] = alpha    # alpha is the largest value in the vector

# %%
# print sparseness of vec_H_DFT
K_DFT, tau_DFT, alpha_DFT = find_real_sparseness(vec_H_AD, threshold=sparse_power_threshold)
sparseness_AD = K_DFT / len(vec_H_AD)
print(
    "K_DFT = {}, \nN_DFT = {}, \n2K_DFT = {} \ntau_DFT = {}, \nalpha_DFT = {}, \nsparseness_DFT = {}".format(K_DFT, len(vec_H_AD), 2*K_DFT, tau_DFT, alpha_DFT, sparseness_AD)
)

# %% [markdown]
# #### Test the NumPy Matrix Indexing

# %%
test_mat = np.array([[1,2,3],[4,5,6],[7,8,9]])
test_mat[0,0]

# %% [markdown]
# #### Test EMRF.py

# %%
from tool_box.EMRF import EMRF
# add noise
noise = np.random.randn(*y.shape) + 1j * np.random.randn(*y.shape)
noise = noise / np.linalg.norm(noise) * np.linalg.norm(y) / 10 ** (SNR / 20) 
y = y + noise

test_EMRF = EMRF(
    Graph = myGraph_AD,
    y = y,
    Phi = Phi_AD,
    vec_H_a=vec_H_AD,
    xi = xi,
    l_x_abs_max=n_x_abs_max,
    l_y_abs_max=n_y_abs_max,
    Geo_Approx_para=Geo_Approx_para,
)


# %%
def test_visual(path = "default"):
    # subplots 1 row 2 columns
    plt.figure(figsize=(14, 5))

    # subplot 1
    plt.subplot(1, 2, 1)
    # visualize the support variable
    visualize_wavenumber_domain_channel(test_EMRF.get_support() + 1)

    # subplot 2
    plt.subplot(1, 2, 2)
    # visualize the hat_h
    visualize_wavenumber_domain_channel(test_EMRF.hat_h_sl)
    # save the figure
    if path != "default":
        plt.savefig(path)


# %% [markdown]
# #### Test $\alpha-\beta$-swap

# %%
iter_num = 20
# use tqdm
for i in tqdm(range(iter_num)):
    test_EMRF.alpha_beta_swap()
    


# %%
NMSE_list_v2 = test_EMRF.NMSE_list_v2
NMSE_list = test_EMRF.NMSE_list
# plot NMSE_list
plt.figure()
plt.subplot(211)
plt.plot(NMSE_list)
# set log scale
plt.yscale('log')
# grid on
plt.grid(True)

plt.subplot(212)
plt.plot(NMSE_list_v2)
# set log scale
plt.yscale('log')
# grid on  
plt.grid(True)

print("NMSE = ", np.min(NMSE_list))
print("NMSE_v2 = ", np.min(NMSE_list_v2))

# %%
# Saving results
import os 
from scipy.io import savemat
# create folder if not exist
# folder = './data', with SNR and N_RF and spacing
spacing = int(lambda_c / delta)
folder = './data' + '/SNR_' + str(SNR) + '_Nx_' + str(N_x) + '_RF_' + str(N_RF) + '_spacing_' + str(spacing)
if not os.path.exists(folder):
    os.makedirs(folder)

test_visual(path=folder + "/RecoveryCompare.png")

savemat(folder + '/alg_GCSE_AD.mat', {
    'NMSE_list_v2' : NMSE_list_v2,
    'NMSE_list' : NMSE_list,
    'meta_data' : {
        # 'cluster_para' : cluster_para,
        'f_c' : f_c,
        'N_x' : N_x,
        'N_y' : N_y,
        'lambda_c' : lambda_c,
        'delta' : delta,
        'N' : N,
        'N_RF' : N_RF,
        'L_x' : L_x,
        'L_y' : L_y,
        'K' : K,
        'sparseness' : sparseness,
        'sparse_power_threshold' : sparse_power_threshold,
        'eta' : eta,
        'xi' : xi,
        'l_x_abs_max' : l_x_abs_max,
        'l_y_abs_max' : l_y_abs_max,
        'M' : M,
        'SNR' : SNR,
        'iter_num' : iter_num
    },
}
)

# %%



