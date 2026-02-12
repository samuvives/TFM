import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

LISTRUNFACTORS = []

x = [1,2,3,4,5]
y = [2,4,8,16,32]

rho, p_valor = spearmanr(x, y)
print(f"ρ = {rho:.3f}, p-valor = {p_valor:.4f}")
