from scipy.optimize import minimize
from prettytable import PrettyTable

# Expected filter efficiencies for when we order filtered DY samples
eff0=0.165
eff1=0.22
eff2=0.21

# Cross-sections in pb for the different jet multiplicities
XS0=1788
XS1=339.7
XS2=125.1
XS_tot = XS0 + XS1 + XS2 


class EarlyRun3:
    # Luminosity in fb-1 (~62.5 fb-1)
    lumi = 7.98 + 26.7 + 18.1 + 9.69
    # Relative uncertainties in final bin
    e0 = 0.4984
    e1 = 0.0376
    e2 = 0.0246
    # Actual content in final bin
    c0 = 0.5033
    c1 = 39.468
    c2 = 54.816
    # Total absolute uncertainty in final bin
    E = (
        ((c0 * e0) ** 2 +  (c1 * e1) ** 2 + (c2 * e2) ** 2) ** 0.5
    )


class Run3_2024:
    # Luminosity in fb-1
    lumi = 109 
    # Relative uncertainties in final bin
    e0 = 0.6234
    e1 = 0.1149
    e2 = 0.1142
    # Actual content in final bin
    c0 = 1.2609
    c1 = 59.906
    c2 = 95.138
    # Total absolute uncertainty in final bin
    E = (
        ((c0 * e0) ** 2 +  (c1 * e1) ** 2 + (c2 * e2) ** 2) ** 0.5
    )
    # Total events originally ordered (NanoAOD event number in millions)
    N_exc_0 = 493.13
    N_exc_1 = 428.39
    N_exc_2 = 232.97
    N_inc = 349.16
    # Total event number broken down by jet multiplicity, assuming the
    # inclusive sample is distributed according to the cross-sections
    N_0 = N_exc_0 + N_inc * (XS0 / XS_tot)
    N_1 = N_exc_1 + N_inc * (XS1 / XS_tot)
    N_2 = N_exc_2 + N_inc * (XS2 / XS_tot)


# target error is EarlyRun3 error scaled by the sqrt of the luminosity ratio
target_E_2024 = EarlyRun3.E * (Run3_2024.lumi / EarlyRun3.lumi) ** 0.5


def constraint(vars):
    x0, x1, x2 = vars
    return (
        ((Run3_2024.e0 * Run3_2024.c0) ** 2 / x0
         + (Run3_2024.e1 * Run3_2024.c1) ** 2 / x1
         + (Run3_2024.e2 * Run3_2024.c2) ** 2 / x2) ** 0.5
         - target_E_2024
    )


def objective(vars):
    x0, x1, x2 = vars
    N0_new = Run3_2024.N_0 * (x0-1) * eff0
    N1_new = Run3_2024.N_1 * (x1-1) * eff1
    N2_new = Run3_2024.N_2 * (x2-1) * eff2
    # we want to minimize the total number of events to be ordered
    return N0_new + N1_new + N2_new


# Initial guesses for x0, x1, x2
initial_guess = [1.0, 1.0, 1.0]

# Define the constraints dictionary
constraints = {
    'type': 'eq',  # Equality constraint
    'fun': constraint
}

# different options for boundaries
bounds_vec = [
    [(1, None), (1, None), (1, None)],
    ]

for bounds in bounds_vec:

    print('\n----------------------------------------------------------------')
    print(f'Performing minimization for bounds: {bounds}')
    # Perform the optimization
    result = minimize(
        objective, initial_guess, bounds=bounds, constraints=constraints
    )

    # Extract results
    x0, x1, x2 = result.x
    print(f"Optimized x0: {x0}, x1: {x1}, x2: {x2}")
    print(f"Objective value: {objective(result.x)}")
    print(f"Constraint value: {constraint(result.x)}")

    # so total number of filtered events to be ordered is the original number
    # of events multiplied by (x-1) and the filter efficiency
    N0_new = Run3_2024.N_0 * (x0-1) * eff0
    N1_new = Run3_2024.N_1 * (x1-1) * eff1
    N2_new = Run3_2024.N_2 * (x2-1) * eff2

    print(
        "\033[91mEvents to be requested for Run3_2024:\033[0m\n"
    )

    table = PrettyTable()
    table.field_names = ["Sample", "Events to be Ordered (millions)"]
    table.add_row(["DYto2Tau_0J", f"{N0_new:.2f}"])
    table.add_row(["DYto2Tau_1J", f"{N1_new:.2f}"])
    table.add_row(["DYto2Tau_2J", f"{N2_new:.2f}"])
    print(table)
