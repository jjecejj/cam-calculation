import customtkinter as ctk
from tkinter import messagebox
from typing import Literal

# Import project classes
try:
    from options import CamSolveOptions, calculate_cam_solve
    from core.cam_geometry.options import CamGeometryOptions
    from vizualization.plotter.options import PlotterOptions
    from vizualization.rotate_animation.options import RotateAnimationOptions
    from exporters.dxf_creator.options import DxfCreatorOptions
except ImportError as e:
    # Fallback for when running directly without setting PYTHONPATH
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from options import CamSolveOptions, calculate_cam_solve
    from core.cam_geometry.options import CamGeometryOptions
    from vizualization.plotter.options import PlotterOptions
    from vizualization.rotate_animation.options import RotateAnimationOptions
    from exporters.dxf_creator.options import DxfCreatorOptions


# Theme settings
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class CamConfiguratorApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Конфигуратор решателя кулачка")
        self.geometry("600x650")

        # Create tabview
        self.tabview = ctk.CTkTabview(self, width=550, height=500)
        self.tabview.pack(padx=20, pady=20, fill="both", expand=True)

        # Add tabs
        self.tab_geom = self.tabview.add("Геометрия")
        self.tab_plot = self.tabview.add("Графики")
        self.tab_anim = self.tabview.add("Анимация")
        self.tab_dxf = self.tabview.add("DXF Экспорт")

        self.init_geometry_tab()
        self.init_plotter_tab()
        self.init_animation_tab()
        self.init_dxf_tab()

        # Generate button
        self.generate_btn = ctk.CTkButton(self, text="Сформировать конфигурацию", command=self.generate_config, height=40)
        self.generate_btn.pack(pady=(0, 20), padx=20, fill="x")

    def init_geometry_tab(self):
        # calculator_type (Literal)
        ctk.CTkLabel(self.tab_geom, text="Тип метода (calculator_type):").grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.geom_calc_type = ctk.CTkOptionMenu(self.tab_geom, values=['polidain', 'polinmail'])
        self.geom_calc_type.grid(row=0, column=1, padx=10, pady=10)

        # kulachok_type (Literal)
        ctk.CTkLabel(self.tab_geom, text="Тип толкателя (kulachok_type):").grid(row=1, column=0, padx=10, pady=10, sticky="w")
        self.geom_kulachok_type = ctk.CTkOptionMenu(self.tab_geom, values=['thin', 'flat', 'roller'])
        self.geom_kulachok_type.grid(row=1, column=1, padx=10, pady=10)

        # N (int)
        ctk.CTkLabel(self.tab_geom, text="Кол-во точек (N):").grid(row=2, column=0, padx=10, pady=10, sticky="w")
        self.geom_n = ctk.CTkEntry(self.tab_geom)
        self.geom_n.insert(0, "1000")
        self.geom_n.grid(row=2, column=1, padx=10, pady=10)

        # initial_angle (float)
        ctk.CTkLabel(self.tab_geom, text="Начальный угол (initial_angle):").grid(row=3, column=0, padx=10, pady=10, sticky="w")
        self.geom_angle = ctk.CTkEntry(self.tab_geom)
        self.geom_angle.insert(0, "0.0")
        self.geom_angle.grid(row=3, column=1, padx=10, pady=10)

        # calculate_optimal_initial_angle (bool)
        self.geom_opt_angle = ctk.CTkCheckBox(self.tab_geom, text="Авторассчет оптимального угла")
        self.geom_opt_angle.select() # Default True
        self.geom_opt_angle.grid(row=4, column=0, columnspan=2, padx=10, pady=10, sticky="w")

    def init_plotter_tab(self):
        # Flags (bool)
        self.plot_tolkatel = ctk.CTkCheckBox(self.tab_plot, text="Графики толкателя")
        self.plot_tolkatel.pack(padx=20, pady=10, anchor="w")

        self.plot_kulachok = ctk.CTkCheckBox(self.tab_plot, text="Графики кулачка")
        self.plot_kulachok.pack(padx=20, pady=10, anchor="w")

        self.plot_profile = ctk.CTkCheckBox(self.tab_plot, text="Показать профиль")
        self.plot_profile.pack(padx=20, pady=10, anchor="w")

        self.plot_together = ctk.CTkCheckBox(self.tab_plot, text="Профиль и графики вместе")
        self.plot_together.pack(padx=20, pady=10, anchor="w")

        # graphs_argument_type (Literal)
        ctk.CTkLabel(self.tab_plot, text="Аргумент графиков:").pack(padx=20, pady=(15, 0), anchor="w")
        self.plot_arg_type = ctk.CTkOptionMenu(self.tab_plot, values=['degree', 'rad', 't'])
        self.plot_arg_type.pack(padx=20, pady=5, anchor="w")

    def init_animation_tab(self):
        self.anim_display = ctk.CTkCheckBox(self.tab_anim, text="Показать анимацию")
        self.anim_display.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        self.anim_together = ctk.CTkCheckBox(self.tab_anim, text="Анимация с графиками")
        self.anim_together.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        self.anim_save = ctk.CTkCheckBox(self.tab_anim, text="Сохранить анимацию")
        self.anim_save.grid(row=2, column=0, padx=10, pady=10, sticky="w")

        self.anim_pause = ctk.CTkCheckBox(self.tab_anim, text="Поддержка паузы")
        self.anim_pause.grid(row=3, column=0, padx=10, pady=10, sticky="w")

        ctk.CTkLabel(self.tab_anim, text="Интервал (мс):").grid(row=0, column=1, padx=10, pady=10, sticky="w")
        self.anim_interval = ctk.CTkEntry(self.tab_anim)
        self.anim_interval.insert(0, "50")
        self.anim_interval.grid(row=0, column=2, padx=10, pady=10)

        ctk.CTkLabel(self.tab_anim, text="Имя файла профиля:").grid(row=1, column=1, padx=10, pady=10, sticky="w")
        self.anim_prof_name = ctk.CTkEntry(self.tab_anim)
        self.anim_prof_name.insert(0, "animation_profile")
        self.anim_prof_name.grid(row=1, column=2, padx=10, pady=10)

        ctk.CTkLabel(self.tab_anim, text="Имя файла дашборда:").grid(row=2, column=1, padx=10, pady=10, sticky="w")
        self.anim_dash_name = ctk.CTkEntry(self.tab_anim)
        self.anim_dash_name.insert(0, "animation_dashboard")
        self.anim_dash_name.grid(row=2, column=2, padx=10, pady=10)

        ctk.CTkLabel(self.tab_anim, text="Аргумент анимации:").grid(row=3, column=1, padx=10, pady=10, sticky="w")
        self.anim_arg_type = ctk.CTkOptionMenu(self.tab_anim, values=['degree', 'rad', 't'])
        self.anim_arg_type.grid(row=3, column=2, padx=10, pady=10)

    def init_dxf_tab(self):
        self.dxf_import = ctk.CTkCheckBox(self.tab_dxf, text="Экспортировать в DXF")
        self.dxf_import.pack(padx=20, pady=20, anchor="w")

        ctk.CTkLabel(self.tab_dxf, text="Имя файла DXF:").pack(padx=20, pady=(5,0), anchor="w")
        self.dxf_name = ctk.CTkEntry(self.tab_dxf, width=200)
        self.dxf_name.insert(0, "kulachok_1")
        self.dxf_name.pack(padx=20, pady=5, anchor="w")

        ctk.CTkLabel(self.tab_dxf, text="Тип геометрии:").pack(padx=20, pady=(15,0), anchor="w")
        self.dxf_line = ctk.CTkOptionMenu(self.tab_dxf, values=["spline", "line"])
        self.dxf_line.pack(padx=20, pady=5, anchor="w")

    def generate_config(self):
        try:
            # Collect data from tabs

            # --- CamGeometryOptions ---
            geom_opts = {
                "calculator_type": self.geom_calc_type.get(),
                "kulachok_type": self.geom_kulachok_type.get(),
                "N": int(self.geom_n.get()),
                "initial_angle": float(self.geom_angle.get()),
                "calculate_optimal_initial_angle": bool(self.geom_opt_angle.get()),
            }

            # --- PlotterOptions ---
            plot_opts = {
                "graphs_tolkatel_flag": bool(self.plot_tolkatel.get()),
                "graphs_kulachok_flag": bool(self.plot_kulachok.get()),
                "graphs_argument_type": self.plot_arg_type.get(),
                "graphs_profile_flag": bool(self.plot_profile.get()),
                "profile_and_graphs_together_flag": bool(self.plot_together.get()),
            }

            # --- RotateAnimationOptions ---
            anim_opts = {
                "display_animation_flag": bool(self.anim_display.get()),
                "animation_profile_and_graphs_together_flag": bool(self.anim_together.get()),
                "save_animation_flag": bool(self.anim_save.get()),
                "animation_intarval": int(self.anim_interval.get()),
                "profile_animation_name_file": self.anim_prof_name.get(),
                "dashboard_animation_name_file": self.anim_dash_name.get(),
                "animation_graphs_argument_type": self.anim_arg_type.get(),
                "animation_pause_flag": bool(self.anim_pause.get()),
            }

            # --- DxfCreatorOptions ---
            dxf_opts = {
                "import_dxf_flag": bool(self.dxf_import.get()),
                "dxf_profile_name": self.dxf_name.get(),
                "dxf_line_type": self.dxf_line.get(),
            }

            # Create instances of Pydantic models
            geom_model = CamGeometryOptions(**geom_opts)
            plot_model = PlotterOptions(**plot_opts)
            anim_model = RotateAnimationOptions(**anim_opts)
            dxf_model = DxfCreatorOptions(**dxf_opts)

            final_config = CamSolveOptions(
                cam_geometry_options=geom_model,
                plotter_options=plot_model,
                rotate_animation_options=anim_model,
                dxf_creator_options=dxf_model
            )

            # Call the solver
            calculate_cam_solve(final_config)

            messagebox.showinfo("Успех", "Расчет успешно выполнен!")

        except ValueError as e:
            messagebox.showerror("Ошибка ввода", f"Пожалуйста, проверьте числовые поля.\nПодробности: {e}")
        except Exception as e:
            messagebox.showerror("Ошибка", str(e))

if __name__ == "__main__":
    app = CamConfiguratorApp()
    app.mainloop()
