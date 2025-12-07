import numpy as np
from tqdm import tqdm
import matplotlib.pylab as plt

def autocorrelation(x):
    x = np.asarray(x)
    x = x - np.mean(x)
    n = len(x)
    f = np.fft.fft(x, n=2*n)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    acf /= acf[0]
    return acf

def pimc_metropolis(N_tau, path, h, m, omega, dtau, hbar=1.0):
    """
    Perform a single PIMC Metropolis sweep over the path.
    Returns updated path and acceptance rate.
    """
    accrate = 0.0
    index = np.random.permutation(N_tau)
    
    # Precompute const ants
    kin_prefactor = m / (2.0 * hbar * dtau)
    pot_prefactor = 0.5 * m * omega**2 * dtau / hbar

    for tau in index:
        # Ensure periodicity
        tau_minus = (tau - 1) % N_tau
        tau_plus = (tau + 1) % N_tau

        # Propose new position
        xnew = path[tau] + h * (np.random.rand() - 0.5)

        # Euclidean action difference
        S_old = kin_prefactor * ((path[tau] - path[tau_minus])**2 + (path[tau_plus] - path[tau])**2) \
                + pot_prefactor * path[tau]**2
        S_new = kin_prefactor * ((xnew - path[tau_minus])**2 + (path[tau_plus] - xnew)**2) \
                + pot_prefactor * xnew**2

        # Metropolis acceptance
        if np.random.rand() < np.exp(-(S_new - S_old)):
            path[tau] = xnew
            accrate += 1 / N_tau

    return path, accrate

def pimc_simulation(N_tau=500, N_sweeps=5000, h=0.05, m=1.0, omega=1.0, 
                    beta=2.0, hbar=1.0, equilibrate=1000, 
                    sample_interval=200):
    """
    Full PIMC simulation for 1D harmonic oscillator with sparse sampling.
    
    Only every `sample_interval` sweeps are recorded to reduce correlations.
    Returns thermal averages <x^2> and <E>.
    """

    dtau = beta / N_tau
    path = np.random.normal(0, 0.1, N_tau)  # initial path
    x2_samples = []
    energy_samples = []
    acceptance_rates = []
    thermal_x2 = []  # store x^2 during thermalization
    
    # ============================================= Thermalization ==========================================
    print("Begin thermalization")
    thermal_x2 = []
    thermal_sample_interval = 10  # record every 10 sweeps for long runs
    tbar = tqdm(range(equilibrate), desc="Thermalization")

    for i in tbar:
        path, acc = pimc_metropolis(N_tau, path, h, m, omega, dtau, hbar)
        
        # Sparse recording
        if i % thermal_sample_interval == 0:
            x2 = np.mean(path**2)
            thermal_x2.append(x2)
            tbar.set_description(f"Step {i+1:05}, acc={acc:.3f}, <x^2>={x2:+8.4f}")

    # Compute running average
    thermal_x2 = np.array(thermal_x2)
    running_avg = np.cumsum(thermal_x2) / np.arange(1, len(thermal_x2)+1)

    # Plot thermalization
    plt.figure(figsize=(8,4))
    plt.plot(thermal_x2, color='lightgray', alpha=0.7, label="Raw ⟨x²⟩")
    plt.plot(running_avg, color='blue', lw=2, label="Running average")
    plt.xlabel("Thermalization step (sampled)")
    plt.ylabel("⟨x²⟩")
    plt.title("Thermalization of ⟨x²⟩ in PIMC")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Measurement sweeps
    print(f"\nMeasuring over {N_sweeps} sweeps (sparse sampling every {sample_interval} sweeps)...")
    mbar = tqdm(range(N_sweeps), desc="Measurement sweeps")
    spare_sweep = 0
    for sweep in mbar:
        path, acc = pimc_metropolis(N_tau, path, h, m, omega, dtau, hbar)
        acceptance_rates.append(acc)
        
        # Only record every sample_interval sweeps
        if sweep % sample_interval == 0:
            spare_sweep += 1
            x2_samples.append(np.mean(path**2))
            kinetic_energy = 0.5 / dtau - (m / (2.0 * dtau**2)) * np.mean((np.roll(path, -1) - path)**2)
            potential_energy = 0.5 * m * omega**2 * np.mean(path**2)
            energy_samples.append(kinetic_energy + potential_energy)
            
            mbar.set_description(f"Sparse Sweep {spare_sweep:05}, <x^2>={x2_samples[-1]:+8.4f}, <E>={energy_samples[-1]:+8.4f}")

    x2_avg = np.mean(x2_samples)
    x2_err = np.std(x2_samples) / np.sqrt(len(x2_samples))
    E_avg = np.mean(energy_samples)
    E_err = np.std(energy_samples) / np.sqrt(len(energy_samples))
    avg_acc_rate = np.mean(acceptance_rates)

    print(f"\nAverage acceptance rate: {avg_acc_rate:.3f}")
    print(f"<x^2> = {x2_avg:.6f} ± {x2_err:.6f}")
    print(f"<E>   = {E_avg:.6f} ± {E_err:.6f}")

    # Compute autocorrelation of x2 samples
    corr = autocorrelation(x2_samples)

    plt.figure(figsize=(7,4))
    plt.plot(corr, lw=2)
    plt.title("Autocorrelation of ⟨x²⟩ samples")
    plt.xlabel("Lag")
    plt.ylabel("Correlation")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return x2_avg, E_avg, x2_err, E_err, avg_acc_rate


# Main block
if __name__ == "__main__":
    N_tau = 100
    N_sweeps = 500000
    h = 0.05                        # initial rate of perturbation
    m = 1.0
    omega = 1.5
    beta = 2.0
    hbar = 1.0 
    thermalization_steps = 50000
    sample_interval = 20
    seed = 0
    np.random.seed(seed)

    print("=" * 50)
    print("PIMC Simulation for 1D Harmonic Oscillator")
    print("=" * 50)
    print(f"Parameters: N_tau={N_tau}, beta={beta}, omega={omega}")

    # Analytical results
    def coth(x): return 1.0 / np.tanh(x)
    analytic_x2 = (hbar / (2.0 * m * omega)) * coth(0.5 * beta * hbar * omega)
    analytic_E = 0.5 * hbar * omega * coth(0.5 * beta * hbar * omega)
    print(f"Analytic <x^2> = {analytic_x2:.6f}")
    print(f"Analytic <E>   = {analytic_E:.6f}")
    print("=" * 50)

    # Adaptive step size tuning (properly decreasing h if acceptance too low)
    print("Tuning step size for ~50% acceptance...")
    test_path = np.random.normal(0, 0.1, N_tau)
    dtau = beta / N_tau
    for _ in range(20):
        test_path, acc_rate = pimc_metropolis(N_tau, test_path, h, m, omega, dtau, hbar)
        if acc_rate > 0.55:
            h *= 1.1
        elif acc_rate < 0.45:
            h *= 0.9
        if 0.45 <= acc_rate <= 0.55:
            break
    print(f"Final step size: {h:.3f}")

    # Run full PIMC
    x2_avg, E_avg, x2_err, E_err, acc_rate = pimc_simulation(
        N_tau=N_tau, N_sweeps=N_sweeps, h=h, m=m, omega=omega, 
        beta=beta, hbar=hbar, equilibrate=thermalization_steps, sample_interval=sample_interval,
    )

    print("\n" + "=" * 50)
    print("RESULTS:")
    print("=" * 50)
    print(f"PIMC <x^2> = {x2_avg:.6f} ± {x2_err:.6f}")
    print(f"Analytic   = {analytic_x2:.6f}")
    print(f"Difference = {abs(x2_avg - analytic_x2):.6f} ({abs(x2_avg - analytic_x2)/analytic_x2*100:.2f}%)")
    print()
    print(f"PIMC <E>   = {E_avg:.6f} ± {E_err:.6f}")
    print(f"Analytic   = {analytic_E:.6f}")
    print(f"Difference = {abs(E_avg - analytic_E):.6f} ({abs(E_avg - analytic_E)/analytic_E*100:.2f}%)")
    print("=" * 50)





