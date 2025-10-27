import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import sqrt, sin, cos, radians

def find_internal_nodes_for_triangle(X, Y, L, angle_deg=120):
    """
    Находит индексы внутренних узлов, лежащих внутри равнобедренного треугольника.
    Треугольник: основание L, боковые стороны = 2, угол при вершине = 120°.
    Возвращает список кортежей (i, j) и словарь node_to_index для нумерации.
    """
    ny, nx = X.shape
    internal_nodes = []
    
    # Параметры треугольника
    angle_rad = radians(angle_deg)
    height = 2 * sin(angle_rad/2)  # Высота треугольника
    
    # Координаты вершин
    A = (0, 0)           # Левая нижняя вершина
    B = (L, 0)           # Правая нижняя вершина  
    C = (L/2, height)    # Верхняя вершина
    
    for j in range(1, ny-1):  # исключаем граничные строки
        for i in range(1, nx-1):  # исключаем граничные столбцы
            x, y = X[j, i], Y[j, i]
            
            # Проверка, находится ли точка внутри треугольника
            if is_point_in_triangle(x, y, A, B, C):
                internal_nodes.append((i, j))
    
    # Создаем словарь: (i, j) -> глобальный индекс
    node_to_index = {node: idx for idx, node in enumerate(internal_nodes)}
    return internal_nodes, node_to_index

def is_point_in_triangle(x, y, A, B, C):
    """
    Проверяет, находится ли точка (x,y) внутри треугольника ABC
    используя метод барицентрических координат.
    """
    x1, y1 = A
    x2, y2 = B  
    x3, y3 = C
    
    # Вычисляем барицентрические координаты
    denom = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    if abs(denom) < 1e-12:
        return False
        
    alpha = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denom
    beta = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denom
    gamma = 1 - alpha - beta
    
    # Точка внутри треугольника, если все координаты между 0 и 1
    return 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1

def build_system_matrix(h, G_theta, X, Y, L, angle_deg=120):
    """
    Строит матрицу A и вектор b для системы A * Phi = b.
    Узлы нумеруются только среди внутренних точек треугольника.
    """
    print(f"\n6. ПОСТРОЕНИЕ СИСТЕМЫ ЛИНЕЙНЫХ УРАВНЕНИЙ (для h={h}):")
    print("   Для внутренних узлов используется шаблон 'крест':")
    print("        φ[i,j] = 1/4 * (φ[i+1,j] + φ[i-1,j] + φ[i,j+1] + φ[i,j-1] + h² * 2 * Gθ)")
    print("   Переносим φ[i,j] в левую часть:")
    print("        -4φ[i,j] + φ[i+1,j] + φ[i-1,j] + φ[i,j+1] + φ[i,j-1] = -h² * 2 * Gθ")
    print("   Или, для удобства:")
    print("        φ[i,j] = 1/4 * (φ[i+1,j] + φ[i-1,j] + φ[i,j+1] + φ[i,j-1]) + (h² * Gθ / 2)")

    # Находим внутренние узлы
    internal_nodes, node_to_index = find_internal_nodes_for_triangle(X, Y, L, angle_deg)
    n = len(internal_nodes)
    
    if n == 0:
        print("   Внутренних узлов не найдено.")
        return None, None, None

    A = np.zeros((n, n))
    b = np.full(n, (h**2 * G_theta) / 2)  # Правая часть: h² * Gθ / 2

    # Заполняем матрицу A
    for idx, (i, j) in enumerate(internal_nodes):
        A[idx, idx] = 1.0  # Диагональный элемент

        # Проверяем соседей: право (i+1, j)
        if (i+1, j) in node_to_index:
            A[idx, node_to_index[(i+1, j)]] = -0.25

        # Лево (i-1, j)
        if (i-1, j) in node_to_index:
            A[idx, node_to_index[(i-1, j)]] = -0.25

        # Верх (i, j+1)
        if (i, j+1) in node_to_index:
            A[idx, node_to_index[(i, j+1)]] = -0.25

        # Низ (i, j-1)
        if (i, j-1) in node_to_index:
            A[idx, node_to_index[(i, j-1)]] = -0.25

    print(f"\n   Матрица коэффициентов A (размер {n}x{n}):")
    for i in range(n):
        row_str = " ".join([f"{A[i, j]:6.3f}" for j in range(n)])
        print(f"   Строка {i}: {row_str} | b[{i}] = {b[i]:.3f}")

    return A, b, node_to_index

def solve_poisson_detailed(h, epsilon, G_theta, L, angle_deg=120):
    """
    Решает уравнение Пуассона в равнобедренном треугольнике.
    """
    # Параметры треугольника
    angle_rad = radians(angle_deg)
    height = 2 * sin(angle_rad/2)  # Высота треугольника
    
    # Ограничивающая прямоугольная область
    Lx, Ly = L, height
    nx = int(Lx / h) + 1
    ny = int(Ly / h) + 1
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    X, Y = np.meshgrid(x, y)

    # Координаты вершин треугольника
    A = (0, 0)           # Левая нижняя вершина
    B = (L, 0)           # Правая нижняя вершина  
    C = (L/2, height)    # Верхняя вершина

    # Инициализация решения
    phi = np.zeros((ny, nx))

    print(f"\n1. ПОСТРОЕНИЕ СЕТКИ:")
    print(f"   Размер ограничивающей области: [{0}, {Lx}] x [{0}, {Ly}]")
    print(f"   Шаг h = {h}")
    print(f"   Количество узлов: nx = {nx}, ny = {ny}")
    print(f"   Параметры треугольника:")
    print(f"     - Основание: L = {L}")
    print(f"     - Боковые стороны: 2")
    print(f"     - Угол при вершине: {angle_deg}°")
    print(f"     - Высота: {height:.3f}")

    # --- Задание граничных условий ---
    print(f"\n2. ЗАДАНИЕ ГРАНИЧНЫХ УСЛОВИЙ φ(x, y) = 0 на границе треугольника:")
    
    # Обнуляем все значения
    phi.fill(0)
    
    # Устанавливаем граничные условия только на границах треугольника
    for j in range(ny):
        for i in range(nx):
            x_val, y_val = X[j, i], Y[j, i]
            
            # Проверка на границу треугольника
            if is_on_triangle_boundary(x_val, y_val, A, B, C):
                phi[j, i] = 0.0

    print(f"   - Все границы треугольника: φ = 0")

    # --- Внутренние узлы: итерационный процесс ---
    print(f"\n3. ИТЕРАЦИОННЫЙ ПРОЦЕСС (метод простой итерации):")
    print(f"   Критерий сходимости: max|φ_new - φ_old| < ε = {epsilon}")
    print(f"   Максимальное число итераций: 10000")
    print(f"   Формула обновления для внутреннего узла (i,j):")
    print(f"        φ[i,j] = 0.25 * (φ[i+1,j] + φ[i-1,j] + φ[i,j+1] + φ[i,j-1] + h² * 2 * Gθ)")
    print(f"   (Применяется только внутри треугольника)")

    max_iter = 10000
    iteration_data = []

    for it in range(max_iter):
        phi_old = phi.copy()
        # Обновление внутренних узлов
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                # Проверка, находится ли узел внутри треугольника
                if is_point_in_triangle(X[j, i], Y[j, i], A, B, C):
                    phi[j, i] = 0.25 * (
                        phi[j, i+1] + phi[j, i-1] +
                        phi[j+1, i] + phi[j-1, i] +
                        h**2 * 2.0 * G_theta
                    )

        iteration_data.append(phi.copy())

        # Проверка сходимости
        diff = np.max(np.abs(phi - phi_old))
        if diff < epsilon:
            print(f"   Сходимость достигнута на итерации {it+1}. Макс. разность = {diff:.6f}")
            break

        # Вывод прогресса каждые 100 итераций
        if (it + 1) % 100 == 0:
            print(f"   Итерация {it+1}: max|Δφ| = {diff:.6f}")

    if it == max_iter - 1:
        print(f"   Предупреждение: Достигнуто максимальное число итераций ({max_iter}). Точность {epsilon} не достигнута.")

    print(f"\n4. РЕЗУЛЬТАТЫ:")
    print(f"   Решение сошлось за {it+1} итераций.")
    
    # Вычисляем min/max только для внутренних точек треугольника
    triangle_phi = []
    for j in range(ny):
        for i in range(nx):
            if is_point_in_triangle(X[j, i], Y[j, i], A, B, C):
                triangle_phi.append(phi[j, i])
    
    if triangle_phi:
        phi_min, phi_max = min(triangle_phi), max(triangle_phi)
        print(f"   Минимальное значение φ внутри треугольника: {phi_min:.6f}")
        print(f"   Максимальное значение φ внутри треугольника: {phi_max:.6f}")
    else:
        print(f"   Внутренних точек треугольника не найдено")

    return X, Y, phi, it+1, iteration_data, y, nx, x, A, B, C

def is_on_triangle_boundary(x, y, A, B, C, tol=1e-10):
    """
    Проверяет, находится ли точка на границе треугольника.
    """
    x1, y1 = A
    x2, y2 = B
    x3, y3 = C
    
    # Проверка на стороне AB
    if abs((y - y1) * (x2 - x1) - (x - x1) * (y2 - y1)) < tol and min(x1, x2) - tol <= x <= max(x1, x2) + tol:
        return True
        
    # Проверка на стороне BC  
    if abs((y - y2) * (x3 - x2) - (x - x2) * (y3 - y2)) < tol and min(x2, x3) - tol <= x <= max(x2, x3) + tol:
        return True
        
    # Проверка на стороне CA
    if abs((y - y3) * (x1 - x3) - (x - x3) * (y1 - y3)) < tol and min(x3, x1) - tol <= x <= max(x3, x1) + tol:
        return True
        
    return False

def create_comparison_table(phi1, phi2, h1, h2, X1, Y1, A, B, C):
    """Создает таблицу сравнения значений φ в общих узлах внутри треугольника."""
    print(f"\n5. СРАВНЕНИЕ РЕШЕНИЙ В ОБЩИХ УЗЛАХ (h={h1} и h={h2}):")
    nx1 = phi1.shape[1]
    ny1 = phi1.shape[0]
    ratio = int(h1 / h2)  # Соотношение шагов

    data = []
    for j in range(ny1):
        for i in range(nx1):
            x_val, y_val = X1[j, i], Y1[j, i]
            # Сравниваем только точки внутри треугольника
            if is_point_in_triangle(x_val, y_val, A, B, C):
                val1 = phi1[j, i]
                val2 = phi2[j * ratio, i * ratio]
                abs_diff = abs(val1 - val2)
                data.append([x_val, y_val, val1, val2, abs_diff])

    if not data:
        print("   Нет общих узлов внутри треугольника для сравнения")
        return None, 0, 0

    df = pd.DataFrame(data, columns=['x', 'y', f'φ (h={h1})', f'φ (h={h2})', 'Абс. Разность'])
    print(df.to_string(index=False, float_format="{:,.4f}".format))

    max_diff = df['Абс. Разность'].max()
    avg_diff = df['Абс. Разность'].mean()
    print(f"\n   Максимальная абсолютная разность: {max_diff:.6f}")
    print(f"   Средняя абсолютная разность: {avg_diff:.6f}")

    return df, max_diff, avg_diff

def plot_triangle(A, B, C, ax, color='blue', linewidth=2):
    """Рисует треугольник на графике."""
    vertices = np.array([A, B, C, A])  # Замыкаем треугольник
    ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=linewidth)

# --- ОСНОВНОЕ ВЫПОЛНЕНИЕ ---

# Параметры из условия
G_theta = 6.0  # Параметр Gθ = 6
L = 2.0        # Основание треугольника
angle_deg = 120 # Угол при вершине
epsilon = 0.1

# Вычисляем высоту треугольника
angle_rad = radians(angle_deg)
height = 2 * sin(angle_rad/2)

print("=" * 60)
print("РЕШЕНИЕ УРАВНЕНИЯ ПУАССОНА В РАВНОБЕДРЕННОМ ТРЕУГОЛЬНИКЕ")
print("=" * 60)
print(f"Параметры задачи:")
print(f"  - Gθ = {G_theta}")
print(f"  - Основание треугольника L = {L}")
print(f"  - Боковые стороны = 2")
print(f"  - Угол при вершине = {angle_deg}°")
print(f"  - Высота треугольника = {height:.3f}")
print(f"  - Точность ε = {epsilon}")

# --- Решение с шагом h = 0.5 ---
X1, Y1, phi1, iter1, iteration_data1, y1, nx1, x1, A, B, C = solve_poisson_detailed(
    0.5, epsilon, G_theta, L, angle_deg)

# --- Построение матрицы системы (для h=0.5) ---
A_matrix, b, node_to_index = build_system_matrix(0.5, G_theta, X1, Y1, L, angle_deg)

# --- Решение с шагом h = 0.25 ---
X2, Y2, phi2, iter2, iteration_data2, y2, nx2, x2, A2, B2, C2 = solve_poisson_detailed(
    0.25, epsilon, G_theta, L, angle_deg)

# --- Сравнение решений ---
comparison_df, max_diff, avg_diff = create_comparison_table(phi1, phi2, 0.5, 0.25, X1, Y1, A, B, C)

# --- ВИЗУАЛИЗАЦИЯ ---
plt.figure(figsize=(18, 14))

# График 1: Область (треугольник) и сетка (h=0.5)
plt.subplot(2, 3, 1)
plot_triangle(A, B, C, plt.gca(), color='blue', linewidth=3)
plt.scatter(X1, Y1, color='gray', s=20, marker='o', label='Узлы сетки', alpha=0.6)

# Закрашиваем внутреннюю область треугольника
triangle_points = np.array([A, B, C])
plt.fill(triangle_points[:, 0], triangle_points[:, 1], 'lightblue', alpha=0.4)

plt.title(f'Равнобедренный треугольник и сетка (h=0.5)\nL={L}, стороны=2, угол={angle_deg}°', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

# График 2: Решение для h = 0.5
plt.subplot(2, 3, 2)
# Создаем маску для треугольника
mask_triangle = np.zeros_like(phi1, dtype=bool)
for j in range(phi1.shape[0]):
    for i in range(phi1.shape[1]):
        mask_triangle[j, i] = is_point_in_triangle(X1[j, i], Y1[j, i], A, B, C)

# Отображаем только значения внутри треугольника
phi1_masked = np.ma.masked_where(~mask_triangle, phi1)
contour1 = plt.contourf(X1, Y1, phi1_masked, levels=20, cmap='viridis')
plt.colorbar(contour1, label='φ(x, y)')
plot_triangle(A, B, C, plt.gca(), color='black', linewidth=2)
plt.title(f'Решение φ(x,y) для h = 0.5\n(Итераций: {iter1})', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, alpha=0.3)

# График 3: Решение для h = 0.25
plt.subplot(2, 3, 3)
mask_triangle2 = np.zeros_like(phi2, dtype=bool)
for j in range(phi2.shape[0]):
    for i in range(phi2.shape[1]):
        mask_triangle2[j, i] = is_point_in_triangle(X2[j, i], Y2[j, i], A, B, C)

phi2_masked = np.ma.masked_where(~mask_triangle2, phi2)
contour2 = plt.contourf(X2, Y2, phi2_masked, levels=20, cmap='viridis')
plt.colorbar(contour2, label='φ(x, y)')
plot_triangle(A, B, C, plt.gca(), color='black', linewidth=2)
plt.title(f'Решение φ(x,y) для h = 0.25\n(Итераций: {iter2})', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, alpha=0.3)

# График 4: Разность решений в общих узлах
plt.subplot(2, 3, 4)
if comparison_df is not None:
    # Создаем сетку для разности
    Z_diff = np.zeros_like(phi1)
    Z_diff.fill(np.nan)  # Заполняем NaN для точек вне треугольника
    
    for j in range(phi1.shape[0]):
        for i in range(phi1.shape[1]):
            if is_point_in_triangle(X1[j, i], Y1[j, i], A, B, C):
                Z_diff[j, i] = abs(phi1[j, i] - phi2[j*2, i*2])
    
    im = plt.contourf(X1, Y1, Z_diff, levels=20, cmap='hot')
    plt.colorbar(im, label='|Δφ|')
    plot_triangle(A, B, C, plt.gca(), color='black', linewidth=2)
    plt.title(f'Абсолютная разность |φ_h=0.5 - φ_h=0.25|\nМакс. = {max_diff:.4f}', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')
else:
    plt.text(0.5, 0.5, 'Нет данных для сравнения', 
             ha='center', va='center', transform=plt.gca().transAxes)

# График 5: Область (треугольник) и сетка (h=0.25)
plt.subplot(2, 3, 5)
plot_triangle(A, B, C, plt.gca(), color='blue', linewidth=3)
plt.scatter(X2, Y2, color='gray', s=10, marker='o', label='Узлы сетки', alpha=0.6)
plt.fill(triangle_points[:, 0], triangle_points[:, 1], 'lightblue', alpha=0.4)
plt.title(f'Равнобедренный треугольник и сетка (h=0.25)\nL={L}, стороны=2, угол={angle_deg}°', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

# График 6: Сходимость (макс. изменение на итерации)
plt.subplot(2, 3, 6)
# Для h=0.5
diffs_h1 = []
phi_prev = iteration_data1[0]
for phi_next in iteration_data1[1:]:
    diff = np.max(np.abs(phi_next - phi_prev))
    diffs_h1.append(diff)
    phi_prev = phi_next

# Для h=0.25
diffs_h2 = []
phi_prev = iteration_data2[0]
for phi_next in iteration_data2[1:]:
    diff = np.max(np.abs(phi_next - phi_prev))
    diffs_h2.append(diff)
    phi_prev = phi_next

plt.semilogy(range(1, len(diffs_h1)+1), diffs_h1, 'b-', label=f'h=0.5 (итераций: {len(diffs_h1)})')
plt.semilogy(range(1, len(diffs_h2)+1), diffs_h2, 'r--', label=f'h=0.25 (итераций: {len(diffs_h2)})')
plt.axhline(y=epsilon, color='green', linestyle='--', label=f'Точность ε={epsilon}')
plt.xlabel('Номер итерации')
plt.ylabel('Макс. |Δφ|')
plt.title('Сходимость метода (логарифмическая шкала)', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# --- ИТОГОВЫЙ АНАЛИЗ И ВЫВОДЫ ---
print(f"\n" + "="*60)
print("ИТОГОВЫЙ АНАЛИЗ")
print("="*60)

print(f"\nКоличество итераций:")
print(f"   - Для h = 0.5: {iter1} итераций")
print(f"   - Для h = 0.25: {iter2} итераций")
print("   Вывод: Уменьшение шага привело к увеличению числа итераций, что ожидаемо из-за роста числа узлов.")

if comparison_df is not None:
    print(f"\nТочность решения:")
    print(f"   - Максимальная абсолютная разность в общих узлах: {max_diff:.6f}")
    print(f"   - Средняя абсолютная разность в общих узлах: {avg_diff:.6f}")
    if max_diff < epsilon:
        print(f"   Вывод: Разность меньше заданной точности ε={epsilon}. Результат с h=0.5 приемлем.")
    else:
        print(f"   Вывод: Разность превышает заданную точность ε={epsilon}. Необходимо использовать меньший шаг (h=0.25 или меньше).")

print(f"\nЭкстремальные значения внутри треугольника:")
# Для h=0.5
triangle_phi1 = []
for j in range(phi1.shape[0]):
    for i in range(phi1.shape[1]):
        if is_point_in_triangle(X1[j, i], Y1[j, i], A, B, C):
            triangle_phi1.append(phi1[j, i])

# Для h=0.25  
triangle_phi2 = []
for j in range(phi2.shape[0]):
    for i in range(phi2.shape[1]):
        if is_point_in_triangle(X2[j, i], Y2[j, i], A, B, C):
            triangle_phi2.append(phi2[j, i])

if triangle_phi1 and triangle_phi2:
    phi1_min, phi1_max = min(triangle_phi1), max(triangle_phi1)
    phi2_min, phi2_max = min(triangle_phi2), max(triangle_phi2)
    print(f"   - Для h=0.5: φ ∈ [{phi1_min:.4f}, {phi1_max:.4f}]")
    print(f"   - Для h=0.25: φ ∈ [{phi2_min:.4f}, {phi2_max:.4f}]")

print(f"\n   Значения функции φ(x,y) в узлах сетки (h=0.5) внутри треугольника:")
print("      x: ", end="")
for i in range(nx1):
    print(f"{x1[i]:5.1f} ", end="")
print()
for j in range(len(y1)-1, -1, -1):
    print(f"y={y1[j]:4.1f}: ", end="")
    for i in range(nx1):
        if is_point_in_triangle(X1[j, i], Y1[j, i], A, B, C):
            print(f"{phi1[j, i]:5.3f} ", end="")
        else:
            print("     - ", end="")
    print()

print(f"\nРЕШЕНИЕ ЗАДАЧИ ДЛЯ РАВНОБЕДРЕННОГО ТРЕУГОЛЬНИКА ЗАВЕРШЕНО")