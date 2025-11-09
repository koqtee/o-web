import numpy as np
import matplotlib.pyplot as plt
from math import radians, sin, cos, sqrt
from scipy.interpolate import griddata
import pandas as pd

def point_in_triangle(x, y, A, B, C, eps=1e-12):
    # barycentric method with numpy
    v0 = np.array([C[0]-A[0], C[1]-A[1]])
    v1 = np.array([B[0]-A[0], B[1]-A[1]])
    v2 = np.array([x-A[0], y-A[1]])
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) < 1e-14:
        return False
    inv = 1.0/denom
    u = (dot11 * dot02 - dot01 * dot12) * inv
    v = (dot00 * dot12 - dot01 * dot02) * inv
    return (u >= -eps) and (v >= -eps) and (u+v <= 1+eps)

def dist_point_to_segment(P, A, B):
    # distance from point P to segment AB
    px, py = P
    ax, ay = A
    bx, by = B
    vx = bx - ax; vy = by - ay
    wx = px - ax; wy = py - ay
    seg_len2 = vx*vx + vy*vy
    if seg_len2 == 0:
        return sqrt((px-ax)**2 + (py-ay)**2)
    t = (vx*wx + vy*wy) / seg_len2
    t = max(0.0, min(1.0, t))
    projx = ax + t*vx
    projy = ay + t*vy
    return sqrt((px-projx)**2 + (py-projy)**2)

def solve_poisson_triangle(h, G_theta=6.0, epsilon=1e-2, omega=1.7, max_iter=20000):
    # Triangle geometry
    L = 2.0
    angle_deg = 120
    angle = radians(angle_deg)
    A = (0.0, 0.0)
    B = (L, 0.0)
    C = (L * cos(angle), L * sin(angle))
    BC_length = sqrt((C[0]-B[0])**2 + (C[1]-B[1])**2)

    # grid covering triangle bbox
    x_min, x_max = min(A[0], B[0], C[0]) - 1e-12, max(A[0], B[0], C[0]) + 1e-12
    y_min, y_max = min(A[1], B[1], C[1]) - 1e-12, max(A[1], B[1], C[1]) + 1e-12
    nx = int(np.floor((x_max - x_min)/h)) + 1
    ny = int(np.floor((y_max - y_min)/h)) + 1
    x = np.linspace(x_min, x_max, nx)
    y = np.linspace(y_min, y_max, ny)
    X, Y = np.meshgrid(x, y)

    # masks
    triangle_mask = np.zeros_like(X, dtype=bool)
    for j in range(ny):
        for i in range(nx):
            triangle_mask[j,i] = point_in_triangle(X[j,i], Y[j,i], A, B, C)

    # boundary identification: nodes in triangle that have any neighbor outside -> boundary nodes
    boundary_node = np.zeros_like(X, dtype=bool)
    for j in range(ny):
        for i in range(nx):
            if not triangle_mask[j,i]: continue
            # check 4-neighbors
            neighbors = [(j, i+1), (j, i-1), (j+1, i), (j-1, i)]
            for jj, ii in neighbors:
                if ii<0 or ii>=nx or jj<0 or jj>=ny or not triangle_mask[jj,ii]:
                    boundary_node[j,i] = True
                    break

    # classify boundary nodes: Dirichlet (on AB) or Neumann (on AC or BC)
    dirichlet_mask = np.zeros_like(X, dtype=bool)
    neumann_mask = np.zeros_like(X, dtype=bool)
    tol = 0.5*h + 1e-12  # tolerance to detect being on segment
    for j in range(ny):
        for i in range(nx):
            if not boundary_node[j,i]: continue
            P = (X[j,i], Y[j,i])
            dAB = dist_point_to_segment(P, A, B)
            dAC = dist_point_to_segment(P, A, C)
            dBC = dist_point_to_segment(P, B, C)
            # Decide: if closest to AB within tol -> Dirichlet; else Neumann.
            if dAB <= dAC + 1e-12 and dAB <= dBC + 1e-12 and dAB < tol:
                dirichlet_mask[j,i] = True
            else:
                neumann_mask[j,i] = True

    # Initialize phi
    phi = np.zeros_like(X)
    # Apply Dirichlet: AB phi=0 (already zero), but keep mask for updates
    phi[dirichlet_mask] = 0.0

    # Iterative SOR scheme with ghost for Neumann (phi_out = phi_in)
    it = 0
    iteration_history = []
    interior_mask = triangle_mask & (~dirichlet_mask)  # nodes to update
    # Precompute neighbor index arrays for speed
    for it in range(1, max_iter+1):
        max_diff = 0.0
        # sweep over interior nodes (including neumann boundary nodes except dirichlet)
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                if not triangle_mask[j,i]: continue  # outside domain
                if dirichlet_mask[j,i]: continue    # fixed value

                # determine neighbor values with handling
                neighbor_vals = []
                ghost_count = 0
                dirichlet_contrib = 0.0

                # right (i+1,j)
                if triangle_mask[j, i+1]:
                    neighbor_vals.append(phi[j, i+1])
                else:
                    # outside: check if it's Neumann boundary (neighboring boundary node is Neumann)
                    # if current node lies on Neumann boundary -> ghost equals current phi (zero normal derivative)
                    # otherwise (if external due to Dirichlet on AB) check neighbor is Dirichlet
                    # Simpler robust approach: find closest point on triangle boundary: if that boundary segment is AB -> dirichlet
                    # We'll use geometric check: point halfway to neighbor lies outside; find which boundary segment is nearest.
                    midx = 0.5*(X[j,i]+X[j,i+1]); midy = 0.5*(Y[j,i]+Y[j,i+1])
                    dAB = dist_point_to_segment((midx,midy), A, B)
                    dAC = dist_point_to_segment((midx,midy), A, C)
                    dBC = dist_point_to_segment((midx,midy), B, C)
                    if dAB <= dAC and dAB <= dBC and dAB < tol:
                        # Dirichlet neighbor value = 0
                        neighbor_vals.append(0.0)
                    else:
                        # Neumann ghost -> phi_out = phi_in
                        neighbor_vals.append(phi[j,i])  # ghost as current
                # left (i-1,j)
                if triangle_mask[j, i-1]:
                    neighbor_vals.append(phi[j, i-1])
                else:
                    midx = 0.5*(X[j,i]+X[j,i-1]); midy = 0.5*(Y[j,i]+Y[j,i-1])
                    dAB = dist_point_to_segment((midx,midy), A, B)
                    dAC = dist_point_to_segment((midx,midy), A, C)
                    dBC = dist_point_to_segment((midx,midy), B, C)
                    if dAB <= dAC and dAB <= dBC and dAB < tol:
                        neighbor_vals.append(0.0)
                    else:
                        neighbor_vals.append(phi[j,i])
                # up (i,j+1)
                if triangle_mask[j+1, i]:
                    neighbor_vals.append(phi[j+1, i])
                else:
                    midx = 0.5*(X[j,i]+X[j+1,i]); midy = 0.5*(Y[j,i]+Y[j+1,i])
                    dAB = dist_point_to_segment((midx,midy), A, B)
                    dAC = dist_point_to_segment((midx,midy), A, C)
                    dBC = dist_point_to_segment((midx,midy), B, C)
                    if dAB <= dAC and dAB <= dBC and dAB < tol:
                        neighbor_vals.append(0.0)
                    else:
                        neighbor_vals.append(phi[j,i])
                # down (i,j-1)
                if triangle_mask[j-1, i]:
                    neighbor_vals.append(phi[j-1, i])
                else:
                    midx = 0.5*(X[j,i]+X[j-1,i]); midy = 0.5*(Y[j,i]+Y[j-1,i])
                    dAB = dist_point_to_segment((midx,midy), A, B)
                    dAC = dist_point_to_segment((midx,midy), A, C)
                    dBC = dist_point_to_segment((midx,midy), B, C)
                    if dAB <= dAC and dAB <= dBC and dAB < tol:
                        neighbor_vals.append(0.0)
                    else:
                        neighbor_vals.append(phi[j,i])

                # Now neighbor_vals has 4 entries. Use standard 5-point Laplacian formula:
                # (phi_E + phi_W + phi_N + phi_S - 4*phi_ij)/h^2 = -G_theta
                # Solve for phi_ij_new:
                sum_neighbors = sum(neighbor_vals)
                rhs = (sum_neighbors + (h**2) * G_theta) / 4.0
                # SOR update: phi_new = (1-omega)*phi_old + omega*rhs
                phi_old = phi[j,i]
                phi_new = (1.0 - omega) * phi_old + omega * rhs
                diff = abs(phi_new - phi_old)
                if diff > max_diff: max_diff = diff
                phi[j,i] = phi_new

        iteration_history.append(max_diff)
        if it % 500 == 0:
            print(f"h={h:.3f}: итерация {it}, maxΔ = {max_diff:.3e}")
        if max_diff < epsilon:
            print(f"h={h:.3f}: сходимость достигнута на итерации {it}, maxΔ = {max_diff:.3e}")
            break

    # post stats
    phi_domain = phi[triangle_mask]
    stats = {
        'iterations': it,
        'phi_min': float(np.min(phi_domain)) if phi_domain.size>0 else 0.0,
        'phi_max': float(np.max(phi_domain)) if phi_domain.size>0 else 0.0,
        'n_nodes': int(np.sum(triangle_mask)),
        'area_mask': float(np.sum(triangle_mask) * (h**2))
    }

    return {
        'X': X, 'Y': Y, 'phi': phi, 'mask': triangle_mask, 'dirichlet': dirichlet_mask,
        'neumann': neumann_mask, 'A': A, 'B': B, 'C': C, 'stats': stats, 'iter_history': iteration_history
    }


G_theta = 6.0
epsilon = 1e-2  # tighter than before
omega = 1.7

res1 = solve_poisson_triangle(0.5, G_theta=G_theta, epsilon=epsilon, omega=omega, max_iter=10000)
res2 = solve_poisson_triangle(0.25, G_theta=G_theta, epsilon=epsilon, omega=omega, max_iter=20000)

# Interpolate coarse -> fine points to compare
points_coarse = np.column_stack([res1['X'][res1['mask']], res1['Y'][res1['mask']]])
values_coarse = res1['phi'][res1['mask']]
points_fine = np.column_stack([res2['X'][res2['mask']], res2['Y'][res2['mask']]])
values_fine = res2['phi'][res2['mask']]

interp_on_fine = griddata(points_coarse, values_coarse, points_fine, method='linear')
# some fine points may have nan due to extrapolation: remove them
valid = ~np.isnan(interp_on_fine)
diffs = np.abs(interp_on_fine[valid] - values_fine[valid])

max_diff = float(np.max(diffs)) if diffs.size>0 else 0.0
avg_diff = float(np.mean(diffs)) if diffs.size>0 else 0.0



# --- ВИЗУАЛИЗАЦИЯ И АНАЛИЗ ---

plt.figure(figsize=(18, 10))

# График 1: Область (треугольник) и сетка (h=0.5)
plt.subplot(2, 3, 1)
X1, Y1, phi1, triangle_mask1 = res1['X'], res1['Y'], res1['phi'], res1['mask']
X2, Y2, phi2, triangle_mask2 = res2['X'], res2['Y'], res2['phi'], res2['mask']
A, B, C = res1['A'], res1['B'], res1['C']
iter1, iter2 = res1['stats']['iterations'], res2['stats']['iterations']
iter_data1, iter_data2 = res1['iter_history'], res2['iter_history']

triangle_area = np.zeros_like(X1, dtype=bool)
for i in range(X1.shape[1]):
    for j in range(X1.shape[0]):
        triangle_area[j, i] = point_in_triangle(X1[j, i], Y1[j, i], A, B, C)

# Отображаем область треугольника
plt.contourf(X1, Y1, triangle_area.astype(int),
             levels=[-0.5, 0.5, 1.5],
             colors=['white', 'lightblue'],
             alpha=0.6)
plt.contour(X1, Y1, triangle_area.astype(int),
            levels=[0.5],
            colors='blue',
            linewidths=2)

# Узлы сетки с шагом 0.5
x_nodes = np.arange(-1.0, 2.1, 0.5)
y_nodes = np.arange(-1.0, 2.1, 0.5)
X_nodes, Y_nodes = np.meshgrid(x_nodes, y_nodes)

# Отображаем только узлы внутри треугольника
for i in range(X_nodes.shape[1]):
    for j in range(X_nodes.shape[0]):
        if point_in_triangle(X_nodes[j, i], Y_nodes[j, i], A, B, C):
            plt.scatter(X_nodes[j, i], Y_nodes[j, i],
                        color='gray', s=25, marker='o', alpha=0.8)

plt.title('Область (треугольник) и сетка (h=0.5)', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.legend(['Граница треугольника', 'Узлы сетки'], loc='upper right')
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.xlim(min(A[0], B[0], C[0]) - 0.1, max(A[0], B[0], C[0]) + 0.1)
plt.ylim(min(A[1], B[1], C[1]) - 0.1, max(A[1], B[1], C[1]) + 0.1)

# График 2: Решение для h=0.5
plt.subplot(2, 3, 2)
phi1_masked = np.where(triangle_mask1, phi1, np.nan)
contour1 = plt.contourf(X1, Y1, phi1_masked, levels=20, cmap='viridis')
plt.colorbar(contour1, label='φ(x, y)')
plt.contour(X1, Y1, phi1_masked, levels=10, colors='black', alpha=0.3)
vertices = np.array([A, B, C, A])
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
plt.title(f'Абсолютная разность решений\nМакс. = {np.nanmax(diff):.4f}', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

# График 5: Сходимость
plt.subplot(2, 3, 5)
plt.semilogy(range(1, len(iter_data1) + 1), iter_data1, 'b-', label=f'h=0.5 ({len(iter_data1)} итерац.)')
plt.semilogy(range(1, len(iter_data2) + 1), iter_data2, 'r--', label=f'h=0.25 ({len(iter_data2)} итерац.)')
plt.axhline(y=epsilon, color='green', linestyle=':', label=f'ε={epsilon}')
plt.xlabel('Номер итерации')
plt.ylabel('max |Δφ|')
plt.title('Сходимость итерационного процесса', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)

# График 6: Область (треугольник) и сетка (h=0.25)
plt.subplot(2, 3, 6)
triangle_area = np.zeros_like(X2, dtype=bool)
for i in range(X2.shape[1]):
    for j in range(X2.shape[0]):
        triangle_area[j, i] = point_in_triangle(X2[j, i], Y2[j, i], A, B, C)

plt.contourf(X2, Y2, triangle_area.astype(int), levels=[-0.5, 0.5, 1.5],
             colors=['white', 'lightblue'], alpha=0.6)
plt.contour(X2, Y2, triangle_area.astype(int), levels=[0.5],
            colors='blue', linewidths=2)

x_nodes = np.arange(-1.0, 2.1, 0.25)
y_nodes = np.arange(-1.0, 2.1, 0.25)
X_nodes, Y_nodes = np.meshgrid(x_nodes, y_nodes)

for i in range(X_nodes.shape[1]):
    for j in range(X_nodes.shape[0]):
        if point_in_triangle(X_nodes[j, i], Y_nodes[j, i], A, B, C):
            plt.scatter(X_nodes[j, i], Y_nodes[j, i],
                        color='gray', s=10, marker='o', alpha=0.7)

plt.title('Область (треугольник) и сетка (h=0.25)', fontsize=12)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.xlim(min(A[0], B[0], C[0]) - 0.1, max(A[0], B[0], C[0]) + 0.1)
plt.ylim(min(A[1], B[1], C[1]) - 0.1, max(A[1], B[1], C[1]) + 0.1)

plt.tight_layout()
plt.show()

# --- ИТОГОВЫЙ АНАЛИЗ ---
max_diff = np.nanmax(diff)
avg_diff = np.nanmean(diff)

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
BC_length = sqrt((C[0] - B[0]) ** 2 + (C[1] - B[1]) ** 2)
print(f"     Стороны: AB = 2, AC = 2, BC = {BC_length:.4f}")
print(f"     Угол A = 120°")
print(f"     Площадь: {0.5 * 2 * 2 * sin(radians(120)):.4f}")

def create_comparison_table_triangle(phi1, phi2, X1, Y1, X2, Y2, domain_mask1, domain_mask2):
    """Сравнение решений в общих узлах для треугольной области"""
    print(f"\n5. СРАВНЕНИЕ РЕШЕНИЙ:")
    
    # Находим общие узлы внутри области (по ближайшим координатам)
    common_points = []
    for i in range(X1.shape[1]):
        for j in range(X1.shape[0]):
            if domain_mask1[j, i]:
                # ищем ближайший узел во второй сетке
                dist = np.sqrt((X2 - X1[j, i])**2 + (Y2 - Y1[j, i])**2)
                min_idx = np.unravel_index(np.argmin(dist), X2.shape)
                
                if domain_mask2[min_idx]:
                    val1 = phi1[j, i]
                    val2 = phi2[min_idx]
                    abs_diff = abs(val1 - val2)
                    common_points.append([X1[j, i], Y1[j, i], val1, val2, abs_diff])
    
    # Создаём таблицу сравнения
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



comparison_df, max_diff, avg_diff = create_comparison_table_triangle(
    phi1, phi2,
    X1, Y1,
    X2, Y2,
    triangle_mask1, triangle_mask2
)
