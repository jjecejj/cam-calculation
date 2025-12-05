import numpy as np
a = np.array([1, 2, 3])

'''
axs[0].plot(t_1, Y_1)
axs[0].plot(t_2, Y_2)
axs[0].plot([t_2[-1], t_3[0]], [r0 + h, r0 + h])
axs[0].plot(t_3, Y_3)
axs[0].plot(t_4, Y_4)
axs[0].set_xlabel('t, c')
axs[0].set_ylabel('координата толкателя, м')
axs[0].grid(True)

# --- 2. Скорость ---
axs[1].plot(t_1, V_1_sek)
axs[1].plot(t_2, V_2_sek)
axs[1].plot([t_2[-1], t_3[0]], [0, 0])
axs[1].plot(t_3, V_3_sek)
axs[1].plot(t_4, V_4_sek)
axs[1].set_xlabel('t, c')
axs[1].set_ylabel('Скорость')
axs[1].grid(True)

# --- 3. Ускорение ---
axs[2].plot(t_1, A_1_sek)
axs[2].plot(t_2, A_2_sek)
axs[2].plot([t_2[-1], t_3[0]], [0, 0])
axs[2].plot(t_3, A_3_sek)
axs[2].plot(t_4, A_4_sek)
axs[2].set_xlabel('t, c')
axs[2].set_ylabel('Ускорение')
axs[2].grid(True)

# --- 4. Рывок ---
axs[3].plot(t_1, D_1_sek)
axs[3].plot(t_2, D_2_sek)
axs[3].plot([t_2[-1], t_3[0]], [0, 0])
axs[3].plot(t_3, D_3_sek)
axs[3].plot(t_4, D_4_sek)
axs[3].set_xlabel('t, c')
axs[3].set_ylabel('Рывок')
axs[3].grid(True)

# --- 5. "Кракен" ---
axs[4].plot(t_1, K_1_sek)
axs[4].plot(t_2, K_2_sek)
axs[4].plot([t_2[-1], t_3[0]], [0, 0])
axs[4].plot(t_3, K_3_sek)
axs[4].plot(t_4, K_4_sek)
axs[4].set_xlabel('t, c')
axs[4].set_ylabel('Кракен')
axs[4].grid(True)

plt.tight_layout()
plt.show()
'''