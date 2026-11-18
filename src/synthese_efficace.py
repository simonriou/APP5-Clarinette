import os
import json
import numpy as np
from scipy.io import wavfile
from numba import njit
from time import time

def charger_parametres(fichier="config.json"):
    with open(fichier, 'r') as f:
        return json.load(f)

def calculer_params_anche(Ks, H, w, r_in, phys):
    rho = phys["rho0"]
    c = phys["c"]
    S = np.pi * (r_in**2)
    Zc = rho * c / S
    pM = Ks * H 
    uA = w * H * np.sqrt(2 * pM / rho) 
    zeta = Zc * uA / pM
    return pM, zeta

@njit
def deriv_F_p_numba(p, gamma, zeta):
    if gamma <= 1e-9:
        return 0.0
    A = zeta * (3.0 * gamma - 1.0) / (2.0 * np.sqrt(gamma))
    B = -zeta * (3.0 * gamma + 1.0) / (8.0 * (gamma ** 1.5))
    C = -zeta * (gamma + 1.0) / (16.0 * (gamma ** 2.5))
    return A + 2.0 * B * p + 3.0 * C * (p ** 2)

@njit
def simuler_boucle_numba_oversampled(omega, alpha, F_n, gamma_array_sim, zeta, fs, duree, oversample):
    N = len(omega)
    dt = 1.0 / (fs * oversample)
    n_steps_out = int(duree * fs)
    n_steps_sim = n_steps_out * oversample
    
    p_out = np.zeros((N, n_steps_out))
    omega2 = omega**2
    p = np.zeros(N)
    v = np.zeros(N)
    
    for i in range(n_steps_sim):
        t = i * dt
        gamma = gamma_array_sim[i]
        
        # RK4
        p_tot1 = np.sum(p)
        v_tot1 = np.sum(v)
        dF1 = deriv_F_p_numba(p_tot1, gamma, zeta)
        dp1 = v
        dv1 = F_n * dF1 * v_tot1 - 2.0 * alpha * v - omega2 * p
        dv1[0] += 1e-2 * np.sin(omega[0] * t)
        
        p2 = p + 0.5 * dt * dp1
        v2 = v + 0.5 * dt * dv1
        dF2 = deriv_F_p_numba(np.sum(p2), gamma, zeta)
        dp2 = v2
        dv2 = F_n * dF2 * np.sum(v2) - 2.0 * alpha * v2 - omega2 * p2
        dv2[0] += 1e-2 * np.sin(omega[0] * (t + 0.5 * dt))
        
        p3 = p + 0.5 * dt * dp2
        v3 = v + 0.5 * dt * dv2
        dF3 = deriv_F_p_numba(np.sum(p3), gamma, zeta)
        dp3 = v3
        dv3 = F_n * dF3 * np.sum(v3) - 2.0 * alpha * v3 - omega2 * p3
        dv3[0] += 1e-2 * np.sin(omega[0] * (t + 0.5 * dt))
        
        p4 = p + dt * dp3
        v4 = v + dt * dv3
        dF4 = deriv_F_p_numba(np.sum(p4), gamma, zeta)
        dp4 = v4
        dv4 = F_n * dF4 * np.sum(v4) - 2.0 * alpha * v4 - omega2 * p4
        dv4[0] += 1e-2 * np.sin(omega[0] * (t + dt))
        
        p = p + (dt / 6.0) * (dp1 + 2.0*dp2 + 2.0*dp3 + dp4)
        v = v + (dt / 6.0) * (dv1 + 2.0*dv2 + 2.0*dv3 + dv4)
        
        if i % oversample == 0:
            p_out[:, i // oversample] = p
            
    return p_out

def simuler_auto_oscillation_depuis_fichier(fichier_freqs, t_env, gamma_env, zeta, fs, duree):
    # Lecture des données modales depuis le fichier CSV
    data = np.loadtxt(fichier_freqs, delimiter=',', skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
        
    N = data.shape[0]
    omega = 2 * np.pi * data[:, 1]
    alpha = data[:, 2]
    F_n = data[:, 4]

    oversample = 4 
    fs_sim = fs * oversample
    
    t_eval_sim = np.arange(0, duree, 1.0 / fs_sim)
    gamma_array_sim = np.interp(t_eval_sim, t_env, gamma_env)
    t_eval_audio = np.arange(0, duree, 1.0 / fs)

    print(f"Intégration RK4 optimisée (suréchantillonnage x{oversample}) sur {duree}s avec {N} modes...")
    
    p_tot_adim = simuler_boucle_numba_oversampled(
        omega, alpha, F_n, gamma_array_sim, zeta, fs, duree, oversample
    )
    
    return t_eval_audio, p_tot_adim

def simulation_et_audio():
    params = charger_parametres("config.json")
    phys = params["physique"]
    fs = params["simulation"]["fs"]
    duree = params["simulation"]["duree"]
    r_in = params["geometrie"]["r_list"][0]
    
    # --- PARAMÉTRAGE DE L'ANCHE ---
    Ks_val = 8e6        
    H_val = 0.4e-3      
    w_val = 12e-3

    start_time = time()
    
    pM, zeta = calculer_params_anche(Ks_val, H_val, w_val, r_in, phys)
    print(f"Paramètres Anche : pM={pM:.1f} Pa, zeta={zeta:.3f}")

    # --- ENVELOPPE DE PRESSION --- # format : (temps en s, pression en Pa)
    enveloppe_souffle_99 = [
        (0.0,  0.0), (0.08, 1200.0), (0.18, 2800.0), (0.35, 2100.0),
        (0.55, 3400.0), (0.75, 2600.0), (0.95, 3200.0), (1.15, 1400.0),
        (1.35, 2400.0), (1.55, 800.0), (1.75, 3000.0), (2.05, 1900.0),
        (2.30, 3600.0), (2.55, 2200.0), (2.75, 2900.0), (3.05, 1000.0),
        (3.35, 2700.0), (3.70, 1700.0), (4.10, 4500.0), (duree, 100.0)       
    ]
    enveloppe_souffle = [
        (0.08, 100.0), (duree / 2, 4000.0), (duree, 0.0)
    ]

    enveloppe_souffle_2 = [
        (0.08, 100.0), (0.1, 2000.0), (1.0, 3500), (duree, 0.0)
    ]
    
    t_env = np.array([pt[0] for pt in enveloppe_souffle])
    pm_env = np.array([pt[1] for pt in enveloppe_souffle])
    gamma_env = pm_env / pM

    # --- SIMULATION ---
    fichier_freqs = "fichier_freqs.txt"
    t, p_modes_adim = simuler_auto_oscillation_depuis_fichier(
        fichier_freqs, t_env, gamma_env, zeta, fs, duree
    )
    
    # Pression interne acoustique dans le bec (Pa)
    p_interne_bec = np.sum(p_modes_adim, axis=0) * pM
    
    # Approximation du rayonnement : dérivation temporelle de la pression du bec
    # p_audio = np.diff(p_interne_bec, prepend=0)
    p_audio = p_interne_bec

    end_time = time()
    print(f"-> Simulation terminée en {end_time - start_time:.2f} secondes.")
    
    # Normalisation et sauvegarde WAV
    max_val = np.max(np.abs(p_audio))
    if max_val > 0:
        p_audio_norm = np.int16((p_audio / max_val) * 32767)
    else:
        p_audio_norm = np.int16(p_audio)
        
    dossier_sortie = "output"
    os.makedirs(dossier_sortie, exist_ok=True)
    fichier_audio = os.path.join(dossier_sortie, "perce_synth_rapide.wav")
    wavfile.write(fichier_audio, fs, p_audio_norm)
    
    print(f"-> Audio généré : {fichier_audio}")

    # # --- AFFICHAGE ---
    # pm_plot = np.interp(t, t_env, pm_env)
    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    
    # ax1.plot(t, pm_plot, color='C1', linewidth=2)
    # ax1.set_title("Profil de Pression de Souffle (Entrée)")
    # ax1.set_ylabel("Pression (Pa)")
    # ax1.grid(True, alpha=0.5)

    # ax2.plot(t, p_interne_bec, color='C4', linewidth=0.8)
    # ax2.set_title("Pression acoustique interne dans le bec (Auto-oscillation)")
    # ax2.set_xlabel("Temps (s)")
    # ax2.set_ylabel("Pression (Pa)")
    # ax2.grid(True, alpha=0.5)
    
    # plt.xlim(0, duree)
    # plt.tight_layout()
    # plt.show()

if __name__ == "__main__":
    simulation_et_audio()