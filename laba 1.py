import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import cos, sin, radians, sqrt

# Функция для проверки принадлежности точки треугольнику (вынесена на верхний уровень)
def point_in_triangle(x, y, A, B, C):
    """Проверяет, лежит ли точка внутри треугольника ABC"""
    # Векторы
    v0 = [C[0] - A[0], C[1] - A[1]]
    v1 = [B[0] - A[0], B[1] - A[1]]
    v2 = [x - A[0], y - A[1]]
    
    # Вычисляем dot products
    dot00 = v0[0]*v0[0] + v0[1]*v0[1]
    dot01 = v0[0]*v1[0] + v0[1]*v1[1]
    dot02 = v0[0]*v2[0] + v0[1]*v2[1]
    dot11 = v1[0]*v1[0] + v1[1]*v1[1]
    dot12 = v1[0]*v2[0] + v1[1]*v2[1]
    
    # Вычисляем барицентрические координаты
    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    
    # Проверяем, находится ли точка внутри треугольника
    return (u >= 0) and (v >= 0) and (u + v <= 1)

def solve_poisson_obtuse_triangle(h, epsilon, G_theta=6):
    """
    Решает уравнение Пуассона в тупоугольном треугольнике 
    с двумя сторонами по 2 и углом 120° между ними
    """
    L = 2.0  # длина двух сторон
    angle_deg = 120  # угол между сторонами длиной 2
    
    # Вершины треугольника:
    # A - вершина с углом 120°
    # AB = AC = 2 (стороны, образующие угол 120°)
    # BC - третья сторона (противоположная углу 120°)
    
    A = (0, 0)  # вершина с углом 120°
    
    # Сторона AB длиной 2 вдоль оси x
    B = (L, 0)
    
    # Сторона AC длиной 2 под углом 120° от AB
    angle_rad = radians(angle_deg)
    C = (L * cos(angle_rad), L * sin(angle_rad))
    
    # Вычисляем длину третьей стороны BC
    BC_length = sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
    
    print(f"Параметры тупоугольного треугольника:")
    print(f"  Стороны: AB = {L}, AC = {L}, угол A = {angle_deg}°")
    print(f"  Третья сторона BC = {BC_length:.4f}")
    print(f"  Вершины: A{A}, B{B}, C{C}")
    
    # Создаем сетку, охватывающую треугольник
    x_min, x_max = min(A[0], B[0], C[0]), max(A[0], B[0], C[0])
    y_min, y_max = min(A[1], B[1], C[1]), max(A[1], B[1], C[1])
    
    nx = int((x_max - x_min) / h) + 1
    ny = int((y_max - y_min) / h) + 1
    x = np.linspace(x_min, x_max, nx)
    y = np.linspace(y_min, y_max, ny)
    X, Y = np.meshgrid(x, y)
    
    # Маска для треугольника
    triangle_mask = np.zeros_like(X, dtype=bool)
    for i in range(nx):
        for j in range(ny):
            triangle_mask[j, i] = point_in_triangle(X[j, i], Y[j, i], A, B, C)
    
    # Маска для границы L (сторона AB длиной 2)
    boundary_mask = np.zeros_like(X, dtype=bool)
    for i in range(nx):
        for j in range(ny):
            if triangle_mask[j, i]:
                # Проверяем, близко ли к стороне AB
                dist_to_AB = abs(Y[j, i])  # AB лежит на оси y=0
                if dist_to_AB < 0.1*h:
                    boundary_mask[j, i] = True
    
    # Инициализация решения
    phi = np.zeros((ny, nx))
    
    print(f"\n1. ПОСТРОЕНИЕ СЕТКИ:")
    print(f"   Область: тупоугольный треугольник")
    print(f"   Стороны: AB = {L}, AC = {L}, угол A = {angle_deg}°")
    print(f"   Третья сторона: BC = {BC_length:.4f}")
    print(f"   Шаг h = {h}")
    print(f"   Количество узлов: nx = {nx}, ny = {ny}")
    print(f"   Узлов в области: {np.sum(triangle_mask)}")
    
    # --- Граничные условия ---
    print(f"\n2. ЗАДАНИЕ ГРАНИЧНЫХ УСЛОВИЙ:")
    print(f"   φ = 0 на стороне L (AB длиной 2)")
    print(f"   ∂φ/∂n = 0 на сторонах AC и BC")
    
    # Устанавливаем граничные условия на стороне L (AB)
    phi[boundary_mask] = 0.0
    
    # --- Итерационный процесс ---
    print(f"\n3. ИТЕРАЦИОННЫЙ ПРОЦЕСС:")
    print(f"   Критерий сходимости: max|φ_new - φ_old| < ε = {epsilon}")
    print(f"   Уравнение: Δφ = -Gθ = -{G_theta}")
    print(f"   Формула обновления для внутренних узлов:")
    print(f"        φ[i,j] = 0.25 * (φ[i+1,j] + φ[i-1,j] + φ[i,j+1] + φ[i,j-1] + h² * Gθ)")
    
    max_iter = 10000
    iteration_data = []
    
    for it in range(max_iter):
        phi_old = phi.copy()
        
        # Обновление только внутренних узлов треугольника
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                if triangle_mask[j, i] and not boundary_mask[j, i]:
                    # Проверяем, что соседи внутри треугольника
                    neighbors = []
                    weights = []
                    
                    # Право (i+1, j)
                    if triangle_mask[j, i+1]:
                        neighbors.append(phi[j, i+1])
                        weights.append(1.0)
                    
                    # Лево (i-1, j)
                    if triangle_mask[j, i-1]:
                        neighbors.append(phi[j, i-1])
                        weights.append(1.0)
                    
                    # Верх (i, j+1)
                    if triangle_mask[j+1, i]:
                        neighbors.append(phi[j+1, i])
                        weights.append(1.0)
                    
                    # Низ (i, j-1)
                    if triangle_mask[j-1, i]:
                        neighbors.append(phi[j-1, i])
                        weights.append(1.0)
                    
                    if len(neighbors) > 0:
                        phi[j, i] = (sum(neighbors) + h**2 * G_theta) / len(neighbors)
        
        iteration_data.append(phi.copy())
        
        # Проверка сходимости
        diff = np.max(np.abs(phi[triangle_mask] - phi_old[triangle_mask]))
        if diff < epsilon:
            print(f"   Сходимость достигнута на итерации {it+1}. Макс. разность = {diff:.6f}")
            break
            
        if (it + 1) % 500 == 0:
            print(f"   Итерация {it+1}: max|Δφ| = {diff:.6f}")
    
    if it == max_iter - 1:
        print(f"   Достигнуто максимальное число итераций ({max_iter})")
    
    # Статистика только по области решения
    phi_in_domain = phi[triangle_mask]
    print(f"\n4. РЕЗУЛЬТАТЫ:")
    print(f"   Итераций: {it+1}")
    print(f"   φ_min = {np.min(phi_in_domain):.6f}")
    print(f"   φ_max = {np.max(phi_in_domain):.6f}")
    print(f"   Узлов в области: {np.sum(triangle_mask)}")
    
    return X, Y, phi, it+1, iteration_data, triangle_mask, boundary_mask, A, B, C

def plot_obtuse_triangle(L=2, angle_deg=120):
    """Визуализация тупоугольного треугольника с двумя сторонами по 2 и углом 120°"""
    angle_rad = radians(angle_deg)
    
    # Вершины треугольника
    A = (0, 0)
    B = (L, 0)
    C = (L * cos(angle_rad), L * sin(angle_rad))
    
    plt.figure(figsize=(10, 8))
    
    # Рисуем треугольник
    vertices = np.array([A, B, C, A])
    plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=2)
    
    # Закрашиваем область
    plt.fill(vertices[:, 0], vertices[:, 1], 'lightblue', alpha=0.5, label='Область решения')
    
    # Подписываем стороны
    plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=3, label='Сторона L (φ=0)')
    plt.plot([A[0], C[0]], [A[1], C[1]], 'g--', linewidth=2, label='Остальные стороны (∂φ/∂n=0)')
    plt.plot([B[0], C[0]], [B[1], C[1]], 'g--', linewidth=2)
    
    # Подписываем вершины
    plt.text(A[0]-0.1, A[1]-0.1, 'A', fontsize=12, ha='right')
    plt.text(B[0]+0.1, B[1]-0.1, 'B', fontsize=12, ha='left')
    plt.text(C[0], C[1]+0.1, 'C', fontsize=12, ha='center')
    
    # Показываем длины сторон и углы
    plt.text(L/2, -0.1, f'L={L}', ha='center', va='top')
    plt.text(L*cos(angle_rad)/3, L*sin(angle_rad)/3, f'L={L}', rotation=angle_deg, ha='center', va='center')
    
    # Показываем угол 120°
    angle_points = 30
    angle_theta = np.linspace(0, angle_rad, angle_points)
    angle_x = 0.3 * np.cos(angle_theta)
    angle_y = 0.3 * np.sin(angle_theta)
    plt.plot(angle_x, angle_y, 'k-', linewidth=1)
    plt.text(0.2*cos(angle_rad/2), 0.2*sin(angle_rad/2), f'{angle_deg}°', ha='center', va='center')
    
    plt.title(f'Тупоугольный треугольник\nAB=AC={L}, угол A={angle_deg}°, Gθ=6', fontsize=14)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    plt.xlim(-0.5, L+0.5)
    plt.ylim(-0.5, L+0.5)
    plt.show()

# --- ОСНОВНОЕ ВЫПОЛНЕНИЕ ---
G_theta = 6
epsilon = 0.001

print("=" * 80)
print("РЕШЕНИЕ УРАВНЕНИЯ ПУАССОНА В ТУПОУГОЛЬНОМ ТРЕУГОЛЬНИКЕ")
print("Уравнение: Δφ = -Gθ")
print(f"Параметры: Gθ = {G_theta}, ε = {epsilon}")
print("Область: тупоугольный треугольник с двумя сторонами по 2 и углом 120° между ними")
print("Граничные условия: φ = 0 на стороне L (AB), ∂φ/∂n = 0 на других сторонах")
print("=" * 80)

# Визуализация области
plot_obtuse_triangle()

# Решение с разными шагами
X1, Y1, phi1, iter1, iter_data1, triangle_mask1, boundary_mask1, A, B, C = solve_poisson_obtuse_triangle(0.5, epsilon, G_theta)
X2, Y2, phi2, iter2, iter_data2, triangle_mask2, boundary_mask2, A, B, C = solve_poisson_obtuse_triangle(0.25, epsilon, G_theta)

def create_comparison_table_triangle(phi1, phi2, X1, Y1, X2, Y2, domain_mask1, domain_mask2):
    """Сравнение решений в общих узлах для треугольной области"""
    print(f"\n5. СРАВНЕНИЕ РЕШЕНИЙ:")
    
    # Находим общие узлы внутри области
    common_points = []
    for i in range(X1.shape[1]):
        for j in range(X1.shape[0]):
            if domain_mask1[j, i]:
                # Находим ближайший узел во второй сетке
                dist = np.sqrt((X2 - X1[j, i])**2 + (Y2 - Y1[j, i])**2)
                min_idx = np.unravel_index(np.argmin(dist), X2.shape)
                
                if domain_mask2[min_idx]:
                    val1 = phi1[j, i]
                    val2 = phi2[min_idx]
                    abs_diff = abs(val1 - val2)
                    common_points.append([X1[j, i], Y1[j, i], val1, val2, abs_diff])
    
    df = pd.DataFrame(common_points, columns=['x', 'y', 'φ (h=0.5)', 'φ (h=0.25)', 'Абс. Разность'])
    
    if len(df) > 0:
        print(f"   Сравнение в {len(df)} общих точках:")
        print(df.head(10).to_string(index=False, float_format="%.4f"))
        if len(df) > 10:
            print("   ... (показаны первые 10 строк)")
        
        max_diff = df['Абс. Разность'].max()
        avg_diff = df['Абс. Разность'].mean()
        print(f"\n   Максимальная разность: {max_diff:.6f}")
        print(f"   Средняя разность: {avg_diff:.6f}")
    else:
        print("   Общие точки для сравнения не найдены")
        max_diff, avg_diff = 0, 0
    
    return df, max_diff, avg_diff

# Сравнение решений
comparison_df, max_diff, avg_diff = create_comparison_table_triangle(phi1, phi2, X1, Y1, X2, Y2, triangle_mask1, triangle_mask2)

# --- ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ ---
plt.figure(figsize=(20, 12))

# Вершины для отрисовки
vertices = np.array([A, B, C, A])

# График 1: Область и сетка (h=0.5)
plt.subplot(2, 3, 1)
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=2)
plt.fill(vertices[:, 0], vertices[:, 1], 'lightblue', alpha=0.3, label='Область решения')

# Подсвечиваем сторону L
plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=3, label='Сторона L (φ=0)')

# Сетка с шагом 0.5 и началом в 0.0
x_nodes = np.arange(-1.0, 2.1, 0.5)  # от 0.0 до 2.0 с шагом 0.5
y_nodes = np.arange(-1.0, 2.1, 0.5)  # от 0.0 до 2.0 с шагом 0.5
X_nodes, Y_nodes = np.meshgrid(x_nodes, y_nodes)

# Отображаем только узлы внутри треугольника
for i in range(7):
    for j in range(7):
        if point_in_triangle(X_nodes[j, i], Y_nodes[j, i], A, B, C):
            plt.scatter(X_nodes[j, i], Y_nodes[j, i], color='gray', s=30, marker='o', alpha=0.7, label='Узлы сетки (h=0.5)' if i == 0 and j == 0 else "")

plt.title('Область решения и сетка (h=0.5)\nТупоугольный треугольник', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

# График 2: Решение для h=0.5
plt.subplot(2, 3, 2)
phi1_masked = np.where(triangle_mask1, phi1, np.nan)
contour1 = plt.contourf(X1, Y1, phi1_masked, levels=20, cmap='viridis')
plt.colorbar(contour1, label='φ(x, y)')
plt.contour(X1, Y1, phi1_masked, levels=10, colors='black', alpha=0.3)
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=1)
plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=2)
plt.title(f'Решение φ(x,y) для h=0.5\n(Итераций: {iter1})', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

# График 3: Решение для h=0.25
plt.subplot(2, 3, 3)
phi2_masked = np.where(triangle_mask2, phi2, np.nan)
contour2 = plt.contourf(X2, Y2, phi2_masked, levels=20, cmap='viridis')
plt.colorbar(contour2, label='φ(x, y)')
plt.contour(X2, Y2, phi2_masked, levels=10, colors='black', alpha=0.3)
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=1)
plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=2)
plt.title(f'Решение φ(x,y) для h=0.25\n(Итераций: {iter2})', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

# График 4: Разность решений
plt.subplot(2, 3, 4)
from scipy.interpolate import griddata
points1 = np.column_stack([X1[triangle_mask1], Y1[triangle_mask1]])
values1 = phi1[triangle_mask1]
points2 = np.column_stack([X2[triangle_mask2], Y2[triangle_mask2]])
phi1_interp = griddata(points1, values1, points2, method='linear')

diff = np.abs(phi1_interp - phi2[triangle_mask2])
diff_plot = np.full_like(phi2, np.nan)
diff_plot[triangle_mask2] = diff

im = plt.contourf(X2, Y2, diff_plot, levels=20, cmap='hot')
plt.colorbar(im, label='|Δφ|')
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=1)
plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=2)
plt.title(f'Абсолютная разность решений\nМакс. = {max_diff:.4f}', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

# График 5: Сходимость
plt.subplot(2, 3, 5)
diffs1 = [np.max(np.abs(iter_data1[i] - iter_data1[i-1])) 
          for i in range(1, len(iter_data1))]
diffs2 = [np.max(np.abs(iter_data2[i] - iter_data2[i-1])) 
          for i in range(1, len(iter_data2))]

plt.semilogy(range(1, len(diffs1)+1), diffs1, 'b-', label=f'h=0.5 ({len(diffs1)} итер.)')
plt.semilogy(range(1, len(diffs2)+1), diffs2, 'r--', label=f'h=0.25 ({len(diffs2)} итер.)')
plt.axhline(y=epsilon, color='green', linestyle=':', label=f'ε={epsilon}')
plt.xlabel('Номер итерации')
plt.ylabel('Макс. |Δφ|')
plt.title('Сходимость итерационного процесса', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)

# График 6: Область и сетка (h=0.25) - аналогично первому графику
plt.subplot(2, 3, 6)
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=2)
plt.fill(vertices[:, 0], vertices[:, 1], 'lightblue', alpha=0.3, label='Область решения')

# Подсвечиваем сторону L
plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=3, label='Сторона L (φ=0)')

# Сетка с шагом 0.25 и началом в 0.0
x_nodes = np.arange(-1.0, 2.1, 0.25)  # от 0.0 до 2.0 с шагом 0.25
y_nodes = np.arange(-1.0, 2.1, 0.25)  # от 0.0 до 2.0 с шагом 0.25
X_nodes, Y_nodes = np.meshgrid(x_nodes, y_nodes)

# Отображаем только узлы внутри треугольника
for i in range(13):
    for j in range(13):
        if point_in_triangle(X_nodes[j, i], Y_nodes[j, i], A, B, C):
            plt.scatter(X_nodes[j, i], Y_nodes[j, i], color='gray', s=20, marker='o', alpha=0.7, label='Узлы сетки (h=0.25)' if i == 0 and j == 0 else "")

plt.title('Область решения и сетка (h=0.25)\nТупоугольный треугольник', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

plt.tight_layout()
plt.show()

# --- ИТОГОВЫЙ АНАЛИЗ ---
print(f"\n6. ИТОГОВЫЙ АНАЛИЗ:")
print(f"   Количество итераций:")
print(f"     h=0.5: {iter1}")
print(f"     h=0.25: {iter2}")

print(f"\n   Точность решения:")
print(f"     Макс. разность: {max_diff:.6f}")
print(f"     Средняя разность: {avg_diff:.6f}")
print(f"     Заданная точность ε: {epsilon}")

if max_diff < epsilon:
    print(f"   ✓ Точность достигнута")
else:
    print(f"   ⚠ Точность не достигнута")

print(f"\n   Экстремальные значения φ в области:")
phi1_domain = phi1[triangle_mask1]
phi2_domain = phi2[triangle_mask2]
print(f"     h=0.5:  min={np.min(phi1_domain):.4f}, max={np.max(phi1_domain):.4f}")
print(f"     h=0.25: min={np.min(phi2_domain):.4f}, max={np.max(phi2_domain):.4f}")

print(f"\n   Параметры треугольника:")
BC_length = sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
print(f"     Стороны: AB = 2, AC = 2, BC = {BC_length:.4f}")
print(f"     Угол A = 120°")
print(f"     Площадь: {0.5 * 2 * 2 * sin(radians(120)):.4f}")
