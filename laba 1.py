import numpy as np
import matplotlib.pyplot as plt

def euler_method(f, x0, y0, b, epsilon):

    h = np.sqrt(epsilon)
    n = int((b - x0) / h) 
    
    x_values = np.linspace(x0, b, n + 1)
    y_values = np.zeros(n + 1)
    y_values[0] = y0
    
    for i in range(n):
        y_values[i + 1] = y_values[i] + h * f(x_values[i], y_values[i])
    
    return x_values, y_values, h

def euler_cauchy_method(f, x0, y0, b, epsilon, max_iter=10):

    h = 0.2
    n = int((b - x0) / h) 
    
    x_values = np.linspace(x0, b, n + 1)
    y_values = np.zeros(n + 1)
    y_values[0] = y0
    
    for i in range(n):
        y_pred = y_values[i] + h * f(x_values[i], y_values[i])
        
        # Итерационное уточнение
        for m in range(1, max_iter + 1):
            y_corrected = y_values[i] + (h / 2) * (
                f(x_values[i], y_values[i]) + 
                f(x_values[i + 1], y_pred)
            )
            
            # Проверка точности
            if abs(y_corrected - y_pred) <= epsilon:
                y_pred = y_corrected
                break
            
            y_pred = y_corrected
        
        y_values[i + 1] = y_pred
    
    return x_values, y_values, h

def get_col(epsilon):
    if epsilon == 0:
        return 6
    col = 0
    while epsilon < 1:
        epsilon *= 10
        col += 1
    return col

def print_results_table(method_name, x_values, y_values, col):
    print(f"\n{method_name}")
    print("k\tx_k\t\ty_k")
    print("-" * 50)

    for k in range(len(x_values)):
        print(f"{k}\t{x_values[k]:.{col}f}\t\t{y_values[k]:.{col}f}")

def solve_differential_equation(f, x0, y0, b, epsilon):

    col = get_col(epsilon)
    
    # 1. Метод Эйлера
    x_euler, y_euler, h_euler = euler_method(f, x0, y0, b, epsilon)
    print(f"Метод Эйлера: шаг h = {h_euler:.1f}, количество шагов = {len(x_euler) - 1}")
    print_results_table("МЕТОД ЭЙЛЕРА:", x_euler, y_euler, col)
    
    # 2. Метод Эйлера-Коши
    x_cauchy, y_cauchy, h_cauchy = euler_cauchy_method(f, x0, y0, b, epsilon)
    print(f"\nМетод Эйлера-Коши: шаг h = {h_cauchy:.1f}, количество шагов = {len(x_cauchy) - 1}")
    print_results_table("МЕТОД ЭЙЛЕРА-КОШИ:", x_cauchy, y_cauchy, col)

def main():
    def f_example1(x, y):
        return 2*x - y
    
    print("y' = 2x - y, y(0) = -1")
    solve_differential_equation(f_example1, x0=0, y0=-1, b=1, epsilon=0.01)
    
    try:
        # Ввод параметров
        func_str = input("Введите функцию f(x, y) для y': ")
        f_user = lambda x, y: eval(func_str)
        x0 = float(input("Введите начальное x0: "))
        y0 = float(input("Введите начальное y0: "))
        b = float(input("Введите конечное b: "))
        epsilon = float(input("Введите требуемую точность ε: "))
        
        # Решение
        solve_differential_equation(f_user, x0, y0, b, epsilon)
        
    except Exception as e:
        print(f"Ошибка: {e}")

if __name__ == "__main__":
    main()