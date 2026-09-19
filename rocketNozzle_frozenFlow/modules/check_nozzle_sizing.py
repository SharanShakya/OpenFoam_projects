# check what nozzle sizing is actually producting
from nozzle_sizing import compute_nozzle_sizing
r = compute_nozzle_sizing(10000, 30e5, 101325, 3.4)
print("Tc:", r["Tc"], "gamma:", r["gamma"], "MW:", r["MW"], "Rt:", r["Rt"], "Re:", r["Re"])