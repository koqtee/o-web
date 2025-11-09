import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import cos, sin, radians, tan, sqrt

def solve_poisson_triangle(h, epsilon, G_theta=6):
    """
    Решает уравнение Пуассона в равнобедренном треугольнике 
    с основанием 2 и углом при вершине 120°
    """
    L = 2.0  # длина основания
    angle_deg = 120  # угол при вершине
    
    # Вычисляем высоту треугольника
    alpha = radians(angle_deg / 2)  # угол между основанием и боковой стороной
    H = (L / 2) * tan(radians(90) - alpha)  # высота треугольника
    
    print(f"Параметры треугольника:")
    print(f"  Основание: {L}")
    print(f"  Угол при вершине: {angle_deg}°")
    print(f"  Высота: {H:.4f}")
    
    # Создаем сетку, охватывающую треугольник
    nx = int(L / h) + 1
    ny = int(H / h) + 1
    x = np.linspace(0, L, nx)
    y = np.linspace(0, H, ny)
    X, Y = np.meshgrid(x, y)
    
    # Маска для треугольника
    # Уравнения сторон треугольника:
    # Левая сторона: y = (H / (L/2)) * x
    # Правая сторона: y = (H / (L/2)) * (L - x)
    triangle_mask = (Y <= (H / (L/2)) * X) & (Y <= (H / (L/2)) * (L - X)) & (Y >= 0)
    
    # Маска для границы L (основание треугольника)
    boundary_mask = (Y < 0.1*h) & triangle_mask  # основание (y ≈ 0)
    
    # Инициализация решения
    phi = np.zeros((ny, nx))
    
    print(f"\n1. ПОСТРОЕНИЕ СЕТКИ:")
    print(f"   Область: равнобедренный треугольник")
    print(f"   Основание: {L}, высота: {H:.4f}, угол: {angle_deg}°")
    print(f"   Шаг h = {h}")
    print(f"   Количество узлов: nx = {nx}, ny = {ny}")
    print(f"   Узлов в области: {np.sum(triangle_mask)}")
    
    # --- Граничные условия ---
    print(f"\n2. ЗАДАНИЕ ГРАНИЧНЫХ УСЛОВИЙ:")
    print(f"   φ = 0 на основании L (y=0)")
    print(f"   ∂φ/∂n = 0 на боковых сторонах")
    
    # Устанавливаем граничные условия на основании
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
        diff = np.max(np.abs(phi - phi_old))
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
    
    return X, Y, phi, it+1, iteration_data, triangle_mask, boundary_mask, L, H

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

def plot_triangle(L=2, angle_deg=120):
    """Визуализация равнобедренного треугольника"""
    # Вычисляем высоту
    alpha = radians(angle_deg / 2)
    H = (L / 2) * tan(radians(90) - alpha)
    
    # Вершины треугольника
    A = (0, 0)  # левая вершина основания
    B = (L, 0)  # правая вершина основания  
    C = (L/2, H)  # вершина
    
    plt.figure(figsize=(8, 6))
    
    # Рисуем треугольник
    vertices = np.array([A, B, C, A])
    plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=2)
    
    # Закрашиваем область
    plt.fill(vertices[:, 0], vertices[:, 1], 'lightblue', alpha=0.5, label='Область решения')
    
    # Подписываем стороны
    plt.plot([A[0], B[0]], [A[1], B[1]], 'r-', linewidth=3, label='Основание (φ=0)')
    plt.plot([A[0], C[0]], [A[1], C[1]], 'g--', linewidth=2, label='Боковая сторона (∂φ/∂n=0)')
    plt.plot([B[0], C[0]], [B[1], C[1]], 'g--', linewidth=2)
    
    plt.title(f'Равнобедренный треугольник\nОснование L={L}, угол={angle_deg}°, Gθ=6', fontsize=14)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    plt.xlim(-0.1, L+0.1)
    plt.ylim(-0.1, H+0.1)
    plt.show()

# --- ОСНОВНОЕ ВЫПОЛНЕНИЕ ---
G_theta = 6
epsilon = 0.001

print("=" * 70)
print("РЕШЕНИЕ УРАВНЕНИЯ ПУАССОНА В РАВНОБЕДРЕННОМ ТРЕУГОЛЬНИКЕ")
print("Уравнение: Δφ = -Gθ")
print(f"Параметры: Gθ = {G_theta}, ε = {epsilon}")
print("Область: равнобедренный треугольник с основанием 2 и углом 120°")
print("Граничные условия: φ = 0 на основании, ∂φ/∂n = 0 на боковых сторонах")
print("=" * 70)

# Визуализация области
plot_triangle()

# Решение с разными шагами
X1, Y1, phi1, iter1, iter_data1, triangle_mask1, boundary_mask1, L, H = solve_poisson_triangle(0.5, epsilon, G_theta)
X2, Y2, phi2, iter2, iter_data2, triangle_mask2, boundary_mask2, L, H = solve_poisson_triangle(0.25, epsilon, G_theta)

# Сравнение решений
comparison_df, max_diff, avg_diff = create_comparison_table_triangle(phi1, phi2, X1, Y1, X2, Y2, triangle_mask1, triangle_mask2)

# --- ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ ---
plt.figure(figsize=(20, 12))

# График 1: Область и сетка (h=0.5)
plt.subplot(2, 3, 1)
# Рисуем треугольник
vertices = np.array([[0, 0], [L, 0], [L/2, H], [0, 0]])
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=2)
plt.fill(vertices[:, 0], vertices[:, 1], 'lightblue', alpha=0.3, label='Область решения')

# Сетка
plt.scatter(X1[triangle_mask1], Y1[triangle_mask1], color='gray', s=20, alpha=0.6, label='Узлы сетки')
plt.title('Область решения и сетка (h=0.5)\nРавнобедренный треугольник', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

# График 2: Решение для h=0.5
plt.subplot(2, 3, 2)
# Создаем маску для отображения только внутри треугольника
phi1_masked = np.where(triangle_mask1, phi1, np.nan)
contour1 = plt.contourf(X1, Y1, phi1_masked, levels=20, cmap='viridis')
plt.colorbar(contour1, label='φ(x, y)')
plt.contour(X1, Y1, phi1_masked, levels=10, colors='black', alpha=0.3)
plt.plot(vertices[:, 0], vertices[:, 1], 'b-', linewidth=1)
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
plt.title(f'Решение φ(x,y) для h=0.25\n(Итераций: {iter2})', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

# График 4: Разность решений
plt.subplot(2, 3, 4)
# Интерполяция phi1 на сетку phi2 для сравнения
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
plt.title(f'Абсолютная разность решений\nМакс. = {max_diff:.4f}', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

# График 5: Сходимость
plt.subplot(2, 3, 5)
# Вычисляем изменения на каждой итерации
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

# График 6: Распределение φ вдоль высоты
plt.subplot(2, 3, 6)
# Выбираем точки вдоль центральной линии (x = L/2)
x_center = L / 2
y_values = np.linspace(0, H, 50)

phi_along_height1 = []
phi_along_height2 = []
for y_val in y_values:
    # Для h=0.5
    dist1 = np.sqrt((X1 - x_center)**2 + (Y1 - y_val)**2)
    idx1 = np.unravel_index(np.argmin(dist1), X1.shape)
    if triangle_mask1[idx1]:
        phi_along_height1.append(phi1[idx1])
    else:
        phi_along_height1.append(np.nan)
    
    # Для h=0.25
    dist2 = np.sqrt((X2 - x_center)**2 + (Y2 - y_val)**2)
    idx2 = np.unravel_index(np.argmin(dist2), X2.shape)
    if triangle_mask2[idx2]:
        phi_along_height2.append(phi2[idx2])
    else:
        phi_along_height2.append(np.nan)

plt.plot(y_values, phi_along_height1, 'bo-', label='h=0.5', markersize=4)
plt.plot(y_values, phi_along_height2, 'r.-', label='h=0.25', markersize=2)
plt.xlabel('Высота y')
plt.ylabel('φ(L/2, y)')
plt.title('Распределение φ вдоль центральной линии', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# --- ВЫВОД ЗНАЧЕНИЙ В УЗЛАХ СЕТКИ ---
print(f"\n6. ЗНАЧЕНИЯ φ(x,y) В УЗЛАХ СЕТКИ (h=0.5):")
print("   (аналогично рис. 6.6 в методичке)")

# Выбираем подмножество узлов для отображения (только в области)
print("\n   Таблица значений (только узлы в треугольнике):")
print("   y\\x", end="")
x_show = np.unique(X1[triangle_mask1])
x_show = x_show[np.argsort(x_show)]
for x_val in x_show[:8]:  # Показываем первые 8 столбцов
    print(f"{x_val:7.1f}", end="")
print()

y_show = np.unique(Y1[triangle_mask1])
y_show = y_show[np.argsort(y_show)]
for y_val in y_show:
    print(f" {y_val:4.1f} ", end="")
    count = 0
    for x_val in x_show:
        if count >= 8:  # Ограничиваем ширину таблицы
            break
        # Находим ближайший узел
        dist = np.sqrt((X1 - x_val)**2 + (Y1 - y_val)**2)
        idx = np.unravel_index(np.argmin(dist), X1.shape)
        if triangle_mask1[idx] and abs(X1[idx] - x_val) < 0.1 and abs(Y1[idx] - y_val) < 0.1:
            print(f"{phi1[idx]:7.3f}", end="")
        else:
            print("       ", end="")
        count += 1
    print()

# --- ИТОГОВЫЙ АНАЛИЗ ---
print(f"\n7. ИТОГОВЫЙ АНАЛИЗ:")
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

print(f"\n   Проверка граничных условий:")
boundary_values = phi1[boundary_mask1]
max_boundary_error = np.max(np.abs(boundary_values))
print(f"     Макс. отклонение на основании (φ=0): {max_boundary_error:.6f}")
if max_boundary_error < 0.01:
    print(f"     ✓ Граничные условия на основании выполнены")
else:
    print(f"     ⚠ Замечание: отклонение на основании")

print(f"\n   Узлов в области:")
print(f"     h=0.5:  {np.sum(triangle_mask1)} узлов")
print(f"     h=0.25: {np.sum(triangle_mask2)} узлов")

print(f"\n   Параметры треугольника:")
print(f"     Основание: {L}")
print(f"     Высота: {H:.4f}")
print(f"     Угол при вершине: 120°")
