import numpy as np
import matplotlib.pyplot as plt

def choose_convenient_step(h_max):
    candidates = [0.1, 0.125, 0.2, 0.25, 0.5]
    for h in reversed(candidates):
        if h <= h_max:
            return h
    return min(candidates)

# Вычисление шага h из неравенства h^5 <= ε
def calculate_step(epsilon):
    h_raw = epsilon ** (1/5)
    return choose_convenient_step(h_raw)

def runge_kutta_4(f, x0, y0, b, h):
    n = int((b - x0) / h) 
    
    x_values = np.linspace(x0, b, n + 1)
    y_values = np.zeros(n + 1)
    y_values[0] = y0
    
        # Списки для хранения коэффициентов
    k1_list = np.zeros(n)
    k2_list = np.zeros(n)
    k3_list = np.zeros(n)
    k4_list = np.zeros(n)

    for i in range(n):
        # Коэффициенты метода Рунге-Кутты
        k1 = h * f(x_values[i], y_values[i])
        k2 = h * f(x_values[i] + h/2, y_values[i] + k1/2)
        k3 = h * f(x_values[i] + h/2, y_values[i] + k2/2)
        k4 = h * f(x_values[i] + h, y_values[i] + k3)
        
        k1_list[i] = k1
        k2_list[i] = k2
        k3_list[i] = k3
        k4_list[i] = k4
        
        # Следующее значение y
        y_values[i + 1] = y_values[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return x_values, y_values, h, k1_list, k2_list, k3_list, k4_list

def get_col(epsilon):
    if epsilon == 0:
        return 6
    col = 0
    while epsilon < 1:
        epsilon *= 10
        col += 1
    return col

def print_results_table(method_name, x_values, y_values, k1_list, k2_list, k3_list, k4_list, col):
    print(f"\n{method_name}")
    print("k\tx_k\t\ty_k\t\tk1\t\tk2\t\tk3\t\tk4")
    print("-" * 95)
    
    # Печатаем начальное условие
    print(f"0\t{x_values[0]:.{col}f}\t\t{y_values[0]:.{col + 1}f}\t\t-\t\t-\t\t-\t\t-")
    
    # Печатаем остальные шаги с коэффициентами
    for k in range(1, len(x_values)):
        print(f"{k}\t{x_values[k]:.{col}f}\t\t{y_values[k]:.{col + 1}f}\t\t" +
              f"{k1_list[k-1]:.{col + 1}f}\t\t{k2_list[k-1]:.{col + 1}f}\t\t" +
              f"{k3_list[k-1]:.{col + 1}f}\t\t{k4_list[k-1]:.{col + 1}f}")
        
def get_user_input():
    
    # Ввод уравнения
    print("Введите дифференциальное уравнение:")
    equation_str = input("y' = ")
    
    # Создаем функцию из строки
    try:
        f = lambda x, y: eval(equation_str)
    except Exception as e:
        print(f"Ошибка в уравнении: {e}")
        return None
    
    # Ввод начальных условий
    x0 = float(input("Начальное значение x0: "))
    y0 = float(input("Начальное значение y0: "))
    b = float(input("Конечное значение x: "))
    epsilon = float(input("Требуемая точность ε: "))
    return f, equation_str, x0, y0, b, epsilon

def solve_differential_equation(f, equation_str, x0, y0, b, h, epsilon):
    
    col = get_col(epsilon)

    # Метод Рунге-Кутты
    x_rk, y_rk, actual_h, k1_list, k2_list, k3_list, k4_list = runge_kutta_4(f, x0, y0, b, h)
    print(f"Шаг: h = {h} Количество шагов: n = {len(x_rk) - 1}")
    
    print_results_table("МЕТОД РУНГЕ-КУТТЫ 4-ГО ПОРЯДКА", x_rk, y_rk, k1_list, k2_list, k3_list, k4_list, col)
    
    return x_rk, y_rk

def main():
    while True:
        # Получаем данные от пользователя
        user_input = get_user_input()
        if user_input is None:
            continue
        
        f, equation_str, x0, y0, b, epsilon = user_input
        
        # Вычисляем шаг h
        h = calculate_step(epsilon)
        
        # Решаем уравнение
        solve_differential_equation(f, equation_str, x0, y0, b, h, epsilon) 

if __name__ == "__main__":
    main()