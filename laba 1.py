import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

def solve_torsion_triangle_by_fem(G_theta, side_length=2.0, angle_deg=120, n_div=8, print_matrices=True):
    """
    Решение задачи кручения треугольной области методом конечных элементов.
    Область: треугольник A(0,0), B(L,0), C(L*cos(angle), L*sin(angle)).
    n_div — число делений стороны (определяет плотность сетки).
    """

    print("="*80)
    print("РЕШЕНИЕ ЗАДАЧИ КРУЧЕНИЯ ТРЕУГОЛЬНИКА МЕТОДОМ КОНЕЧНЫХ ЭЛЕМЕНТОВ")
    print(f"Параметр Gθ = {G_theta}")
    print(f"Сторона L = {side_length}, угол между сторонами = {angle_deg}°")
    print("="*80)

    # --- 1. ГЕНЕРАЦИЯ УЗЛОВ СЕТКИ ---
    print("\n1. ГЕНЕРАЦИЯ УЗЛОВ СЕТКИ:")

    L = side_length
    angle = np.radians(angle_deg)
    A = np.array([0.0, 0.0])
    B = np.array([L, 0.0])
    C = np.array([L * np.cos(angle), L * np.sin(angle)])

    n = n_div
    nodes = []
    idx_map = {}

    for i in range(n + 1):
        for j in range(n + 1 - i):
            pt = (i * B + j * C) / n
            idx = len(nodes)
            nodes.append(pt)
            idx_map[(i, j)] = idx
    nodes = np.array(nodes)
    n_nodes = len(nodes)

    print(f"   Узлов всего: {n_nodes}")
    if print_matrices:
        print("   № узла |     x     |     y    ")
        print("   -------|-----------|-----------")
        for i, (x, y) in enumerate(nodes):
            print(f"   {i+1:6d} | {x:9.3f} | {y:9.3f}")

    # --- 2. ТРИАНГУЛЯЦИЯ ---
    print("\n2. ТРИАНГУЛЯЦИЯ ОБЛАСТИ:")

    elements = []
    for i in range(n):
        for j in range(n - i):
            n1 = idx_map[(i, j)]
            n2 = idx_map[(i + 1, j)]
            n3 = idx_map[(i, j + 1)]
            elements.append([n1, n2, n3])
            if j + i + 1 < n:
                n4 = idx_map[(i + 1, j + 1)]
                elements.append([n2, n4, n3])
    elements = np.array(elements)
    n_elems = len(elements)

    print(f"   Элементов создано: {n_elems}")

    # --- 3. МАТРИЦЫ ЖЁСТКОСТИ ---
    print("\n3. ПОСТРОЕНИЕ МАТРИЦ ЖЁСТКОСТИ:")

    def triangle_area(p1, p2, p3):
        return 0.5 * abs((p2[0]-p1[0])*(p3[1]-p1[1]) - (p3[0]-p1[0])*(p2[1]-p1[1]))

    def element_stiffness_matrix(elem_nodes):
        p1, p2, p3 = elem_nodes
        area = triangle_area(p1, p2, p3)
        xi, yi = p1
        xj, yj = p2
        xk, yk = p3
        bi, bj, bk = yj - yk, yk - yi, yi - yj
        ci, cj, ck = xk - xj, xi - xk, xj - xi
        B = np.array([[bi, bj, bk], [ci, cj, ck]]) / (2 * area)
        K_e = area * (B.T @ B)
        return K_e, area

    K_global = np.zeros((n_nodes, n_nodes))
    areas = []

    for e_idx, elem in enumerate(elements):
        pts = nodes[elem]
        K_e, area = element_stiffness_matrix(pts)
        areas.append(area)
        for i_local, i_global in enumerate(elem):
            for j_local, j_global in enumerate(elem):
                K_global[i_global, j_global] += K_e[i_local, j_local]
        if print_matrices and e_idx < 3:
            print(f"\n   Элемент {e_idx+1}, узлы {elem}: площадь = {area:.4f}")
            print("   Локальная матрица K_e:")
            for r in range(3):
                print("      " + " ".join(f"{K_e[r,c]:10.6f}" for c in range(3)))

    # --- 4. ВЕКТОР НАГРУЗКИ ---
    print("\n4. ПОСТРОЕНИЕ ВЕКТОРА НАГРУЗКИ:")

    F_global = np.zeros(n_nodes)
    for e_idx, elem in enumerate(elements):
        f_e = (G_theta * areas[e_idx] / 3) * np.ones(3)
        for i_local, i_global in enumerate(elem):
            F_global[i_global] += f_e[i_local]

    if print_matrices:
        print("   Первые 10 значений F:")
        for i in range(min(10, n_nodes)):
            print(f"   F[{i+1}] = {F_global[i]:.6f}")

    # --- 5. ГРАНИЧНЫЕ УСЛОВИЯ ---
    print("\n5. ГРАНИЧНЫЕ УСЛОВИЯ:")
    boundary_nodes = [i for i, (x, y) in enumerate(nodes) if abs(y) < 1e-8]
    print(f"   Граничных узлов (φ=0): {len(boundary_nodes)}")

    K_original = K_global.copy()
    F_original = F_global.copy()

    for node in boundary_nodes:
        K_global[node, :] = 0.0
        K_global[:, node] = 0.0
        K_global[node, node] = 1.0
        F_global[node] = 0.0

    if print_matrices:
        print("   Модифицированная матрица K (первые 5×5):")
        for i in range(min(5, n_nodes)):
            print("   " + " ".join(f"{K_global[i,j]:8.3f}" for j in range(min(5, n_nodes))))

    # --- 6. РЕШЕНИЕ СИСТЕМЫ ---
    print("\n6. РЕШЕНИЕ СИСТЕМЫ:")

    try:
        Phi = np.linalg.solve(K_global, F_global)
        print("   Система успешно решена")
    except np.linalg.LinAlgError:
        print("   Ошибка: матрица вырождена")
        return None

    if print_matrices:
        print("\n   Первые 10 значений φ:")
        for i in range(min(10, n_nodes)):
            print(f"   φ[{i+1}] = {Phi[i]:.6f}")

    # --- 7. КРУТЯЩИЙ МОМЕНТ ---
    print("\n7. КРУТЯЩИЙ МОМЕНТ:")

    integral_phi = 0.0
    for e_idx, elem in enumerate(elements):
        avg_phi = np.mean(Phi[elem])
        integral_phi += areas[e_idx] * avg_phi
    T = 2 * G_theta * integral_phi

    print(f"   Крутящий момент T = {T:.6f}")

    # --- 8. ВИЗУАЛИЗАЦИЯ ---
    print("\n8. ВИЗУАЛИЗАЦИЯ:")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    tri = Triangulation(nodes[:, 0], nodes[:, 1], elements)

    # Сетка
    for elem in elements:
        pts = nodes[elem]
        ax1.plot([pts[0,0], pts[1,0], pts[2,0], pts[0,0]],
                 [pts[0,1], pts[1,1], pts[2,1], pts[0,1]], 'b-', lw=1)
    ax1.scatter(nodes[:,0], nodes[:,1], color='r', s=30)
    ax1.set_title("Сетка конечных элементов треугольника")
    ax1.axis('equal')
    ax1.grid(True)

    # Поле φ
    cf = ax2.tricontourf(tri, Phi, levels=20, cmap='viridis')
    ax2.tricontour(tri, Phi, levels=10, colors='k', linewidths=0.4)
    plt.colorbar(cf, ax=ax2, label='φ(x,y)')
    ax2.set_title("Распределение функции φ(x,y)")
    ax2.axis('equal')
    ax2.grid(True)

    # Матрица K до ГУ
    im3 = ax3.imshow(np.abs(K_original), cmap='hot', aspect='auto')
    plt.colorbar(im3, ax=ax3)
    ax3.set_title("Матрица жесткости K (до ГУ)")

    # Матрица K после ГУ
    im4 = ax4.imshow(np.abs(K_global), cmap='hot', aspect='auto')
    plt.colorbar(im4, ax=ax4)
    ax4.set_title("Матрица жесткости K (после ГУ)")

    plt.tight_layout()
    plt.show()

    # --- 9. ИТОГОВАЯ ТАБЛИЦА ---
    print("\n" + "="*80)
    print("ИТОГОВЫЕ РЕЗУЛЬТАТЫ")
    print("="*80)
    print(f"Параметр Gθ = {G_theta}")
    print(f"Сторона треугольника L = {side_length}")
    print(f"Угол = {angle_deg}°")
    print(f"Крутящий момент T = {T:.6f}")
    print(f"Количество узлов = {n_nodes}")
    print(f"Количество элементов = {n_elems}")

    print("\nЗНАЧЕНИЯ ФУНКЦИИ φ(x,y) В УЗЛАХ:")
    print("№ узла |     x     |     y     |    φ(x,y)    ")
    print("-------|-----------|-----------|--------------")
    for i in range(n_nodes):
        print(f"{i+1:6d} | {nodes[i,0]:9.3f} | {nodes[i,1]:9.3f} | {Phi[i]:12.6f}")

    return nodes, elements, Phi, T, K_original, F_original, K_global, F_global


# --- ПРИМЕР ЗАПУСКА ---
if __name__ == "__main__":
    G_theta = 6.0
    nodes, elements, Phi, T, K_original, F_original, K_bc, F_bc = \
        solve_torsion_triangle_by_fem(G_theta, side_length=2.0, angle_deg=120, n_div=8, print_matrices=True)
