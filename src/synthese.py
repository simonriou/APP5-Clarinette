import os
import json
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import scipy.signal
from scipy.io import wavfile
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def charger_parametres(fichier="config.json"):
    with open(fichier, 'r') as f:
        return json.load(f)

def p0_func(x, params_p0):
    return params_p0["amplitude"] * np.exp(-0.5 * ((x - params_p0["x0"]) / params_p0["sigma"])**2)

def dp0_dt_func(x, params_dp0):
    return params_dp0["amplitude"]

def calcul_Zc_kc_eff(omega, r1, r2, params_phys, params_dissip):
    rho0 = params_phys["rho0"]
    c = params_phys["c"]
    gamma = params_phys["gamma"]
    Cp = params_phys["Cp"]
    
    rm = (r1 + r2) / 2.0
    if not params_dissip["activer"]:
        return rho0 * c, omega / c

    if params_dissip["mode_entree"] == "coefficients":
        mu = params_dissip["mu"]
        kappa = params_dissip["kappa"]
    else:
        omega_ref = 2 * np.pi * params_dissip["f_ref"]
        delta_v = params_dissip["delta_v"]
        delta_t = params_dissip["delta_t"]
        mu = (delta_v**2) * rho0 * omega_ref
        kappa = (delta_t**2) * rho0 * Cp * omega_ref

    lv = mu / (rho0 * c)
    Pr = (mu * Cp) / kappa
    
    rv = rm * np.sqrt(omega / (c * lv))
    rt = rv * np.sqrt(Pr)
    sqrt_mj = np.sqrt(-1j)
    
    rho_s = rho0 * (1.0 + 2.0 * sqrt_mj / rv)
    chi_s = 1.0 / (rho0 * c**2)
    chi_t = chi_s * (1.0 + 2.0 * (gamma - 1.0) * sqrt_mj / rt)
    
    k = omega * np.sqrt(rho_s * chi_t)
    rho_c_eff = np.sqrt(rho_s / chi_t)
    return rho_c_eff, k

def matrice_transfert_segment(omega, L, r1, r2, phys, dissip):
    rho_c_eff, k = calcul_Zc_kc_eff(omega, r1, r2, phys, dissip)
    
    if np.abs(r1 - r2) < 1e-9:
        Zc = rho_c_eff / (np.pi * r1**2)
        return np.array([
            [np.cos(k * L), 1j * Zc * np.sin(k * L)],
            [1j * (1.0 / Zc) * np.sin(k * L), np.cos(k * L)]
        ], dtype=complex)
    
    x1 = r1 * L / (r2 - r1)
    x2 = r2 * L / (r2 - r1)
    
    M11 = (r2 / r1) * np.cos(k * L) - np.sin(k * L) / (k * x1)
    M12 = 1j * rho_c_eff / (np.pi * r1 * r2) * np.sin(k * L)
    M21 = (np.pi * r1 * r2 / rho_c_eff) * (
        1j * np.sin(k * L) * (1.0 + 1.0 / (k**2 * x1 * x2)) 
        + (np.cos(k * L) / (1j * k)) * (1.0 / x1 - 1.0 / x2)
    )
    M22 = (r1 / r2) * np.cos(k * L) + np.sin(k * L) / (k * x2)
    
    return np.array([[M11, M12], [M21, M22]], dtype=complex)

def calcul_Zr(omega, r_out, phys, ray):
    if ray["type"] == "aucun":
        return 0j
        
    rho0 = phys["rho0"]
    c = phys["c"]
    S = np.pi * r_out**2
    k0 = omega / c
    negliger_Rr = ray["negliger_Rr"]

    if ray["type"] == "baffle":
        Rr = ((k0 * r_out)**2) / 2.0 if not negliger_Rr else 0.0
        Xr = (8.0 * k0 * r_out) / (3.0 * np.pi)
    else:
        Rr = ((k0 * r_out)**2) / 4.0 if not negliger_Rr else 0.0
        Xr = 0.6133 * k0 * r_out
        
    return (rho0 * c / S) * (Rr + 1j * Xr)

def calcul_facteur_modal_Fn(phi_entree, r_entree, c, N_n):
    if N_n == 0:
        return 0.0
    S_entree = np.pi * r_entree**2
    F_n = c * S_entree * (phi_entree ** 2) / N_n
    return F_n

def calculer_params_anche(Ks, H, w, r_in, phys):
    rho = phys["rho0"]
    c = phys["c"]
    S = np.pi * (r_in**2)
    Zc = rho * c / S
    pM = Ks * H # pression de plaquage
    uA = w * H * np.sqrt(2 * pM / rho) # paramètre de débit
    zeta = Zc * uA / pM
    
    return pM, zeta

def deriv_F_p(p, gamma, zeta):
    """ Modèle polynomial (ordre 3) pour une stabilité numérique garantie """
    # Sécurité physique et numérique : pas de souffle = pas de couplage
    if gamma <= 1e-9:
        return 0.0
        
    # Calcul des coefficients
    A = zeta * (3.0 * gamma - 1.0) / (2.0 * np.sqrt(gamma))
    B = -zeta * (3.0 * gamma + 1.0) / (8.0 * (gamma ** 1.5))
    C = -zeta * (gamma + 1.0) / (16.0 * (gamma ** 2.5))
    
    # Dérivée polynomiale : F'(p) = A + 2Bp + 3Cp^2
    return A + 2.0 * B * p + 3.0 * C * (p ** 2)

def simuler_auto_oscillation(resultats_modaux, t_env, gamma_env, zeta, fs, duree):
    N = len(resultats_modaux)

    omega = np.zeros(N)
    alpha = np.zeros(N)
    F_n = np.zeros(N)

    for i, res in enumerate(resultats_modaux):
        omega[i] = 2 * np.pi * res[1]
        alpha[i] = res[2]
        F_n[i] = res[4]

    def systeme_dyn(t, y):
        p = y[:N]
        v = y[N:]
        
        p_tot = np.sum(p)
        v_tot = np.sum(v)
        
        # Interpolation linéaire dynamique de la pression de souffle au temps t
        gamma_actuel = np.interp(t, t_env, gamma_env)
            
        deriv_F = deriv_F_p(p_tot, gamma_actuel, zeta)
        
        dp_dt = v
        dv_dt = F_n * deriv_F * v_tot - 2 * alpha * v - (omega**2) * p
        
        # Simulation de la turbulence constante du souffle
        bruit_turbulence = 1e-2 * np.sin(omega[0] * t)
        dv_dt[0] += bruit_turbulence 
        
        return np.concatenate((dp_dt, dv_dt))

    y0 = np.zeros(2 * N)
    t_eval = np.arange(0, duree, 1.0 / fs)

    print(f"Intégration ODE sur {duree}s avec {N} modes. . .")

    sol = integrate.solve_ivp(
        systeme_dyn,
        [0, duree],
        y0,
        t_eval=t_eval,
        method='Radau'
    )

    if not sol.success:
        print(f"ATTENTION : Le solveur a échoué à t={sol.t[-1]:.3f}s. Motif : {sol.message}")

    p_tot_adim = sol.y[:N, :]
    return sol.t, p_tot_adim

def simulation_et_audio():

    def pression_forcee(z_eval, f_vals, z_list, r_list, phys, dissip, ray):
        n_seg = len(z_list) - 1
        r_out = r_list[-1]
        seg = np.clip(np.searchsorted(z_list, z_eval, side="right") - 1, 0, n_seg - 1)
        r_eval = r_list[seg] + (z_eval - z_list[seg]) / (z_list[seg + 1] - z_list[seg]) * (r_list[seg + 1] - r_list[seg])

        P = np.zeros((len(f_vals), len(z_eval)), dtype=complex)

        for a, f in enumerate(f_vals):
            omega = 2 * np.pi * f
            M_seg = [matrice_transfert_segment(omega, z_list[i+1] - z_list[i],
                                            r_list[i], r_list[i+1], phys, dissip)
                    for i in range(n_seg)]
            M_tot = np.eye(2, dtype=complex)
            for M in M_seg:
                M_tot = M_tot @ M

            Zr = calcul_Zr(omega, r_out, phys, ray)
            s = 1.0 / (M_tot[1, 0] * Zr + M_tot[1, 1])
            X_node = [np.array([Zr, 1.0], dtype=complex) * s]
            for i in range(n_seg - 1, -1, -1):
                X_node.insert(0, M_seg[i] @ X_node[0])

            for b, z in enumerate(z_eval):
                i = seg[b]
                if np.isclose(z, z_list[i + 1]):
                    P[a, b] = X_node[i + 1][0]
                else:
                    M_part = matrice_transfert_segment(omega, z_list[i + 1] - z, r_eval[b],
                                                    r_list[i + 1], phys, dissip)
                    P[a, b] = (M_part @ X_node[i + 1])[0]
        return P

    def plot_pression_3d(z_eval, f_vals, P, echelle="lineaire", elev0=30, azim0=-60):
        amp = np.abs(P)
        if echelle == "dB":
            amp = 20 * np.log10(np.maximum(amp, 1e-12))
            zlabel = "|P| (dB, réf. 1 Pa·s·m⁻³)"
        else:
            zlabel = "|P| (Pa·s·m⁻³)"

        z_mesh, f_mesh = np.meshgrid(z_eval, f_vals)
        fig = plt.figure(figsize=(11, 8))
        ax = fig.add_subplot(111, projection="3d")
        surf = ax.plot_surface(z_mesh, f_mesh, amp, cmap="viridis",
                            linewidth=0, antialiased=True,
                            rcount=min(400, len(f_vals)), ccount=len(z_eval))
        fig.colorbar(surf, ax=ax, pad=0.1, shrink=0.7, label=zlabel)
        ax.set_title("Pression dans la perce (débit d'entrée unitaire)")
        ax.set_xlabel("Coordonnée spatiale z (m)")
        ax.set_ylabel("Fréquence (Hz)")
        ax.set_zlabel(zlabel)
        ax.view_init(elev=elev0, azim=azim0)

        fig.subplots_adjust(bottom=0.17)
        ax_azim = fig.add_axes([0.20, 0.08, 0.60, 0.025])
        ax_elev = fig.add_axes([0.20, 0.03, 0.60, 0.025])
        s_azim = Slider(ax_azim, "Azimut (°)", -180, 180, valinit=azim0)
        s_elev = Slider(ax_elev, "Élévation (°)", -90, 90, valinit=elev0)

        def maj_vue(_):
            ax.view_init(elev=s_elev.val, azim=s_azim.val)
            fig.canvas.draw_idle()

        s_azim.on_changed(maj_vue)
        s_elev.on_changed(maj_vue)

        def sync_curseurs(event):
            azim = (ax.azim + 180) % 360 - 180
            elev = max(-90, min(90, ax.elev))
            for s, val in ((s_azim, azim), (s_elev, elev)):
                s.eventson = False
                s.set_val(val)
                s.eventson = True
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("button_release_event", sync_curseurs)
        fig._curseurs = (s_azim, s_elev)
        plt.show()
        return fig

    params = charger_parametres("config.json")
    phys = params["physique"]
    c = phys["c"]
    
    z_list = np.array(params["geometrie"]["z_list"])
    r_list = np.array(params["geometrie"]["r_list"])
    L_tot = z_list[-1] - z_list[0]
    r_out = r_list[-1]
    
    n_modes = params["simulation"]["n_modes"]
    fs = params["simulation"]["fs"]
    duree = params["simulation"]["duree"]
    ray = params["rayonnement"]
    dissip = params["dissipation"]
    
    dossier_sortie = "output"
    os.makedirs(dossier_sortie, exist_ok=True)
    nom_base = f"perce_auto_oscillation_{len(z_list)-1}segments"
    fichier_audio = os.path.join(dossier_sortie, f"{nom_base}.wav")
    fichier_freqs = os.path.join(dossier_sortie, f"{nom_base}_resultats.txt")

    resultats_modaux = []
    z_plot = np.linspace(z_list[0], z_list[-1], 300)
    
    print(f"Balayage fréquentiel pour trouver les pics sur {len(z_list)-1} segments...")
    f_scan = np.arange(20.0, 2500.0, 0.5) 
    mag_Zin_scan = np.zeros_like(f_scan)
    
    for idx, f in enumerate(f_scan):
        omega = 2 * np.pi * f
        M_tot = np.eye(2, dtype=complex)
        for i in range(len(z_list) - 1):
            M_tot = M_tot @ matrice_transfert_segment(omega, z_list[i+1]-z_list[i], r_list[i], r_list[i+1], phys, dissip)
        Zr = calcul_Zr(omega, r_out, phys, ray)
        Zin_val = (M_tot[0, 0] * Zr + M_tot[0, 1]) / (M_tot[1, 0] * Zr + M_tot[1, 1])
        mag_Zin_scan[idx] = np.abs(Zin_val)
        
    mag_Zin_scan = np.nan_to_num(mag_Zin_scan, nan=0.0)

    pics_idx, _ = scipy.signal.find_peaks(mag_Zin_scan)
    pics_idx = pics_idx[:n_modes] 
    
    if len(pics_idx) == 0:
        print("Erreur : Aucun pic de résonance trouvé.")
        return

    print(f"{len(pics_idx)} modes identifiés. Calcul des racines complexes et facteurs modaux...")
    
    for n, idx_pic in enumerate(pics_idx, 1):
        f_guess = f_scan[idx_pic]
        omega_0 = 2 * np.pi * f_guess
        
        def systeme_ferme(omega):
            M_tot = np.eye(2, dtype=complex)
            for i in range(len(z_list) - 1):
                M_tot = M_tot @ matrice_transfert_segment(omega, z_list[i+1]-z_list[i], r_list[i], r_list[i+1], phys, dissip)
            Zr = calcul_Zr(omega, r_out, phys, ray)
            return M_tot[1, 0] * Zr + M_tot[1, 1]
            
        omega_complexe = optimize.newton(systeme_ferme, omega_0 + 0j, tol=1e-6, maxiter=200, disp=False)
        omega_rn = np.real(omega_complexe)
        sigma_n = np.abs(np.imag(omega_complexe))
        
        M_tot_reel = np.eye(2, dtype=complex)
        for i in range(len(z_list) - 1):
            M_tot_reel = M_tot_reel @ matrice_transfert_segment(omega_rn, z_list[i+1]-z_list[i], r_list[i], r_list[i+1], phys, dissip)
        Zr_reel = calcul_Zr(omega_rn, r_out, phys, ray)
        Zin_val = (M_tot_reel[0, 0] * Zr_reel + M_tot_reel[0, 1]) / (M_tot_reel[1, 0] * Zr_reel + M_tot_reel[1, 1])
        mag_Zin = np.abs(Zin_val)
        
        X_frontieres = [np.array([calcul_Zr(omega_complexe, r_out, phys, ray), 1.0], dtype=complex)]
        
        for i in range(len(z_list) - 1, 0, -1):
            L_i = z_list[i] - z_list[i-1]
            M_i = matrice_transfert_segment(omega_complexe, L_i, r_list[i-1], r_list[i], phys, dissip)
            X_frontieres.insert(0, M_i @ X_frontieres[0])

        max_P = max(np.abs(X[0]) for X in X_frontieres)
        if max_P > 0:
            X_frontieres = [X / max_P for X in X_frontieres]
            
        def r_func(z):
            i = min(np.searchsorted(z_list, z, side='right') - 1, len(z_list) - 2)
            i = max(0, i)
            return r_list[i] + (z - z_list[i]) / (z_list[i+1] - z_list[i]) * (r_list[i+1] - r_list[i])
            
        def eval_phi_n_complex(z):
            i = min(np.searchsorted(z_list, z, side='right') - 1, len(z_list) - 2)
            i = max(0, i)
            if np.isclose(z, z_list[i+1], rtol=0.0, atol=1e-12): return X_frontieres[i+1][0]
            M_part = matrice_transfert_segment(omega_complexe, z_list[i+1] - z, r_func(z), r_list[i+1], phys, dissip)
            return (M_part @ X_frontieres[i+1])[0]

        z_grid = np.linspace(0, L_tot, 50)
        phi_grid = [eval_phi_n_complex(z_val) for z_val in z_grid]
        phase_correction = np.angle(phi_grid[np.argmax(np.abs(phi_grid))])
            
        f_max = f_scan[-1]
        n_pts = max(2000, int(40 * L_tot * f_max / c))
        z_fine = np.linspace(z_list[0], z_list[-1], n_pts)
        phi_c  = np.array([eval_phi_n_complex(z) for z in z_fine])
        phi    = np.real(phi_c * np.exp(-1j * phase_correction))
        S_fine = np.pi * np.array([r_func(z) for z in z_fine])**2

        N_n = integrate.simpson(S_fine * phi**2, x=z_fine)
            
        if N_n == 0: continue

        # Pression modale réelle au bec (z = 0)
        phi_0 = np.real(eval_phi_n_complex(0.0) * np.exp(-1j * phase_correction))

        # Dérivée spatiale en sortie, via le débit modal : dp/dz = -j*omega*rho0*U/S
        # (X_frontieres est déjà normalisé par max_P, comme phi)
        S_out = np.pi * r_out**2
        U_L = X_frontieres[-1][1]
        dphi_dx_L = np.real(
            -1j * omega_complexe * phys["rho0"] * U_L / S_out * np.exp(-1j * phase_correction)
        )

        F_n = calcul_facteur_modal_Fn(phi_0, r_list[0], c, N_n)
        resultats_modaux.append((n, omega_rn / (2 * np.pi), sigma_n, mag_Zin, F_n, phi_0, dphi_dx_L))

    # --- PARAMÉTRAGE DE L'ANCHE ---
    Ks_val = 8e6        # Raideur surfacique (Pa/m)
    H_val = 0.4e-3      # Ouverture (m)
    w_val = 12e-3       # Largeur du canal (m)
    
    # Calcul des constantes physiques de l'anche (pM et zeta ne dépendent pas de pm)
    pM, zeta = calculer_params_anche(Ks_val, H_val, w_val, r_list[0], phys)
    print(f"Paramètres Anche : pM={pM:.1f} Pa, zeta={zeta:.3f}")

    # --- ENVELOPPE DE PRESSION PERSONNALISÉE ---
    # Points de la courbe sous la forme (Temps en secondes, Pression en Pa)
    enveloppe_souffle_test = [
        (0.0,   0.0),      # Départ à 0
        (0.05,  2500.0),   # Attaque rapide : montée à 2500 Pa
        (0.5,   2500.0),   # Maintien
        (0.8,   1800.0),   # Decrescendo vers 1800 Pa
        (1.2,   1800.0),   # Maintien
        (1.8,   100.0),    # Baisse
        (2.5,   100.0),    # Maintien
        (3.0,   3500.0),   # Remontée
        (duree, 0.0)       # Garantie de couverture jusqu'à la fin
    ]

    enveloppe_souffle = [
        (0.0,  0.0),       # Entrée progressive
        (0.08, 1200.0),
        (0.18, 2800.0),
        (0.35, 2100.0),
        (0.55, 3400.0),
        (0.75, 2600.0),
        (0.95, 3200.0),
        (1.15, 1400.0),
        (1.35, 2400.0),
        (1.55, 800.0),
        (1.75, 3000.0),
        (2.05, 1900.0),
        (2.30, 3600.0),
        (2.55, 2200.0),
        (2.75, 2900.0),
        (3.05, 1000.0),
        (3.35, 2700.0),
        (3.70, 1700.0),
        (4.10, 3200.0),
        (duree, 0.0)       # Relâchement en fin de note
]
    
    # Extraction et adimensionnement de l'enveloppe
    t_env = np.array([pt[0] for pt in enveloppe_souffle])
    pm_env = np.array([pt[1] for pt in enveloppe_souffle])
    gamma_env = pm_env / pM

    # --- SIMULATION NON LINÉAIRE ---
    t, p_modes_adim = simuler_auto_oscillation(
        resultats_modaux, 
        t_env, 
        gamma_env, 
        zeta, 
        fs, 
        duree
    )
    
    # Paramètres de rayonnement et d'écoute
    D = 2.0 if ray["type"] == "baffle" else 4.0
    r_obs = 1.0 # Distance d'écoute en mètres (modifiable)
    facteur_rad = (np.pi * r_out**2) / (D * np.pi * r_obs)
    delai = r_obs / c
    
    p_audio_physique = np.zeros_like(t)
    p_interne_bec = np.zeros_like(t) # On garde la pression bec pour l'affichage (optionnel)
    
    for i, res in enumerate(resultats_modaux):
        phi_0 = res[5]
        dphi_dx_L = res[6]
        
        # Pression interne physique du mode au niveau du bec (Pa)
        p_mode_bec = p_modes_adim[i] * pM
        p_interne_bec += p_mode_bec
        
        # Amplitude modale temporelle absolue q_n(t)
        q_n_t = p_mode_bec / phi_0
        
        # Application du retard lié à la distance de propagation
        t_retard = t - delai
        q_n_t_retard = np.interp(t_retard, t, q_n_t, left=0.0)
        
        # Sommation de la contribution rayonnée du mode i
        p_audio_physique += -facteur_rad * dphi_dx_L * q_n_t_retard

    # Normalisation et sauvegarde WAV (inchangé)
    max_val = np.max(np.abs(p_audio_physique))
    if max_val > 0:
        p_audio_norm = np.int16((p_audio_physique / max_val) * 32767)
    else:
        p_audio_norm = np.int16(p_audio_physique)
        
    wavfile.write(fichier_audio, fs, p_audio_norm)
    
    # Sauvegarde des résultats modaux
    with open(fichier_freqs, 'w') as f:
        # Vous pouvez choisir d'ajouter phi_0 et dphi_dx_L dans l'en-tête si vous voulez les conserver
        f.write("Mode, Frequence(Hz), Taux_Amort(s^-1), Module_Zin(Pa.s.m^-3), Facteur_Modal_Fn\n")
        
        for n, f_n, sigma_n, mag_Zin, F_n, phi_0, dphi_dx_L in resultats_modaux:
            f.write(f"{n}, {f_n:.3f}, {sigma_n:.3f}, {mag_Zin:.3f}, {F_n:.5e}\n")

    print(f"-> Audio généré : {fichier_audio}")
    print(f"-> Résultats sauvegardés : {fichier_freqs}")

    # --- AFFICHAGE DES RÉSULTATS TEMPORELS ---
    # Interpolation de l'enveloppe pour l'affichage avec la même base de temps que l'audio
    pm_plot = np.interp(t, t_env, pm_env)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    
    ax1.plot(t, pm_plot, color='C1', linewidth=2)
    ax1.set_title("Profil de Pression de Souffle (Entrée)")
    ax1.set_ylabel("Pression dans la bouche (Pa)")
    ax1.grid(True, alpha=0.5)

    ax2.plot(t, p_audio_physique, color='C4', linewidth=0.8)
    ax2.set_title("Pression acoustique dans le bec (Auto-oscillation)")
    ax2.set_xlabel("Temps (s)")
    ax2.set_ylabel("Pression (Pa)")
    ax2.grid(True, alpha=0.5)
    
    plt.xlim(0, duree)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    simulation_et_audio()