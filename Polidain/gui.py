import json
import customtkinter as ctk
from tkinter import messagebox, filedialog
from typing import Dict, Any

try:
    from options import CamSolveOptions, calculate_cam_solve
    from core.cam_geometry.options import CamGeometryOptions
    from core.cam_geometry.config import KulachokConfig
    from core.profiling_methods.polidain.config import PolidainConfig
    from core.profiling_methods.polinmail.config import PolinmailConfig, LocalPolinmailConfig
    from vizualization.plotter.options import PlotterOptions
    from vizualization.rotate_animation.options import RotateAnimationOptions
    from exporters.dxf_creator.options import DxfCreatorOptions
except ImportError:
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from options import CamSolveOptions, calculate_cam_solve
    from core.cam_geometry.options import CamGeometryOptions
    from core.cam_geometry.config import KulachokConfig
    from core.profiling_methods.polidain.config import PolidainConfig
    from core.profiling_methods.polinmail.config import PolinmailConfig, LocalPolinmailConfig
    from vizualization.plotter.options import PlotterOptions
    from vizualization.rotate_animation.options import RotateAnimationOptions
    from exporters.dxf_creator.options import DxfCreatorOptions

ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")


class CamConfiguratorApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Конфигуратор решателя кулачка")
        self.geometry("750x850")

        # Виджет с вкладками
        self.tabview = ctk.CTkTabview(self, width=700, height=700)
        self.tabview.pack(padx=20, pady=10, fill="both", expand=True)

        self.tab_geom = self.tabview.add("Общие")
        self.tab_kulachok = self.tabview.add("Кулачок")
        self.tab_method = self.tabview.add("Метод расчета")
        self.tab_plot = self.tabview.add("Графики")
        self.tab_anim = self.tabview.add("Анимация")
        self.tab_dxf = self.tabview.add("DXF Экспорт")

        self.entries = {}

        self.init_geometry_tab()
        self.init_kulachok_tab()
        self.init_method_tab()
        self.init_plotter_tab()
        self.init_animation_tab()
        self.init_dxf_tab()

        # --- ИНИЦИАЛИЗАЦИЯ СОСТОЯНИЙ ---
        self.toggle_initial_angle()
        self.update_follower_params(self.geom_kulachok_type.get())
        self.toggle_dxf_name()

        # --- Панель кнопок внизу ---
        self.btn_frame = ctk.CTkFrame(self)
        self.btn_frame.pack(pady=(0, 20), padx=20, fill="x")

        self.btn_load = ctk.CTkButton(self.btn_frame, text="Загрузить JSON", command=self.load_from_json,
                                      fg_color="gray50")
        self.btn_load.pack(side="left", padx=10, expand=True, fill="x")

        self.btn_save = ctk.CTkButton(self.btn_frame, text="Сохранить JSON", command=self.save_to_json,
                                      fg_color="gray50")
        self.btn_save.pack(side="left", padx=10, expand=True, fill="x")

        self.btn_gen = ctk.CTkButton(self.btn_frame, text="Сформировать конфиг", command=self.generate_config)
        self.btn_gen.pack(side="left", padx=10, expand=True, fill="x")

    def create_entry(self, parent, key, label_text, default_val, row, col=0):
        ctk.CTkLabel(parent, text=label_text).grid(row=row, column=col * 2, padx=10, pady=5, sticky="w")
        entry = ctk.CTkEntry(parent, width=150)
        entry.insert(0, str(default_val))
        entry.grid(row=row, column=col * 2 + 1, padx=10, pady=5, sticky="w")
        self.entries[key] = entry

    # --- ФУНКЦИИ БЛОКИРОВКИ ПОЛЕЙ ---
    def toggle_initial_angle(self):
        """Блокирует поле 'Начальный угол', если включен авторассчет"""
        if self.geom_opt_angle.get() == 1:
            self.entries["geom_angle"].configure(state="disabled")
        else:
            self.entries["geom_angle"].configure(state="normal")

    def update_follower_params(self, choice):
        """Блокирует параметры толкателя в зависимости от его типа"""
        if choice == 'thin':
            self.entries["k_Dt"].configure(state="disabled")
            self.entries["k_Rr"].configure(state="disabled")
        elif choice == 'flat':
            self.entries["k_Dt"].configure(state="normal")
            self.entries["k_Rr"].configure(state="disabled")
        elif choice == 'roller':
            self.entries["k_Dt"].configure(state="normal")
            self.entries["k_Rr"].configure(state="normal")

    def toggle_dxf_name(self):
        """Блокирует поле имени DXF, если экспорт выключен"""
        if self.dxf_import.get() == 1:
            self.entries["dxf_name"].configure(state="normal")
        else:
            self.entries["dxf_name"].configure(state="disabled")

    # --- ИНИЦИАЛИЗАЦИЯ ВКЛАДОК ---
    def init_geometry_tab(self):
        ctk.CTkLabel(self.tab_geom, text="Тип метода:").grid(row=0, column=0, padx=10, pady=5, sticky="w")
        self.geom_calc_type = ctk.CTkOptionMenu(self.tab_geom, values=['polidain', 'polinmail'],
                                                command=self.update_method_visibility)
        self.geom_calc_type.grid(row=0, column=1, padx=10, pady=5)

        ctk.CTkLabel(self.tab_geom, text="Тип толкателя:").grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.geom_kulachok_type = ctk.CTkOptionMenu(self.tab_geom, values=['thin', 'flat', 'roller'],
                                                    command=self.update_follower_params)
        self.geom_kulachok_type.grid(row=1, column=1, padx=10, pady=5)

        self.create_entry(self.tab_geom, "geom_n", "Кол-во точек (N):", 1000, 2)
        self.create_entry(self.tab_geom, "geom_angle", "Начальный угол:", 0.0, 3)

        self.geom_opt_angle = ctk.CTkCheckBox(self.tab_geom, text="Авторассчет оптимального угла",
                                              command=self.toggle_initial_angle)
        self.geom_opt_angle.select()
        self.geom_opt_angle.grid(row=4, column=0, columnspan=2, padx=10, pady=10, sticky="w")

    def init_kulachok_tab(self):
        self.create_entry(self.tab_kulachok, "k_Nk", "Обороты в минуту (N_k):", 100, 0)
        self.create_entry(self.tab_kulachok, "k_D", "Базовый диаметр (D):", 0.05, 1)
        self.create_entry(self.tab_kulachok, "k_h", "Макс. перемещение (h):", 0.02, 2)
        self.create_entry(self.tab_kulachok, "k_z", "Тепловой зазор (z):", 0.001, 3)
        self.create_entry(self.tab_kulachok, "k_Dt", "Диаметр толкателя (D_t):", 0.0, 4)
        self.create_entry(self.tab_kulachok, "k_Rr", "Радиус ролика (R_r):", 0.0, 5)

        self.create_entry(self.tab_kulachok, "k_fpod", "Фаза подъёма (рад):", 1.0, 0, col=1)
        self.create_entry(self.tab_kulachok, "k_fv", "Фаза выдержки (рад):", 1.0, 1, col=1)
        self.create_entry(self.tab_kulachok, "k_fop", "Фаза опускания (рад):", 1.0, 2, col=1)
        self.create_entry(self.tab_kulachok, "k_fz", "Фаза зазора (рад):", 0.1, 3, col=1)

    def init_method_tab(self):
        self.frame_polidain = ctk.CTkFrame(self.tab_method, fg_color="transparent")
        self.frame_polinmail = ctk.CTkFrame(self.tab_method, fg_color="transparent")

        # Polidain Fields
        ctk.CTkLabel(self.frame_polidain, text="Настройки Polidain", font=("Arial", 16, "bold")).grid(row=0, column=0,
                                                                                                      columnspan=2,
                                                                                                      pady=(10, 5),
                                                                                                      sticky="w")
        self.create_entry(self.frame_polidain, "pd_m", "Степень m (>=2):", 3, 1)
        self.create_entry(self.frame_polidain, "pd_d", "Разность d (>=1):", 1, 2)
        self.create_entry(self.frame_polidain, "pd_k1", "k_1 (агресс. зазор):", 2, 3)
        self.create_entry(self.frame_polidain, "pd_k2", "k_2 (агресс. подъем):", 2, 4)
        self.create_entry(self.frame_polidain, "pd_k3", "k_3 (агресс. опуск.):", 2, 5)
        self.create_entry(self.frame_polidain, "pd_k4", "k_4 (агресс. зазор):", 2, 6)

        # Polinmail Fields
        ctk.CTkLabel(self.frame_polinmail, text="Настройки Polinmail (Config 1)", font=("Arial", 16, "bold")).grid(
            row=0, column=0, columnspan=2, pady=(10, 5), sticky="w")
        self.create_entry(self.frame_polinmail, "pm_m", "Степень m (>=1):", 1, 1)
        self.create_entry(self.frame_polinmail, "pm_d", "Разность d (>=1):", 1, 2)
        self.create_entry(self.frame_polinmail, "pm_bc", "Граничные условия (через запятую):", "-1, 0, 0, 0, 0", 3)

        self.frame_polidain.pack(fill="both", expand=True)

    def update_method_visibility(self, choice):
        if choice == "polidain":
            self.frame_polinmail.pack_forget()
            self.frame_polidain.pack(fill="both", expand=True)
        elif choice == "polinmail":
            self.frame_polidain.pack_forget()
            self.frame_polinmail.pack(fill="both", expand=True)

    def init_plotter_tab(self):
        self.plot_tolkatel = ctk.CTkCheckBox(self.tab_plot, text="Графики толкателя")
        self.plot_tolkatel.pack(padx=20, pady=10, anchor="w")

        self.plot_kulachok = ctk.CTkCheckBox(self.tab_plot, text="Графики кулачка")
        self.plot_kulachok.pack(padx=20, pady=10, anchor="w")

        self.plot_profile = ctk.CTkCheckBox(self.tab_plot, text="Показать профиль")
        self.plot_profile.pack(padx=20, pady=10, anchor="w")

        self.plot_together = ctk.CTkCheckBox(self.tab_plot, text="Профиль и графики вместе")
        self.plot_together.pack(padx=20, pady=10, anchor="w")

        ctk.CTkLabel(self.tab_plot, text="Аргумент графиков:").pack(padx=20, pady=(15, 0), anchor="w")
        self.plot_arg_type = ctk.CTkOptionMenu(self.tab_plot, values=['degree', 'rad', 't'])
        self.plot_arg_type.pack(padx=20, pady=5, anchor="w")

    def init_animation_tab(self):
        self.anim_display = ctk.CTkCheckBox(self.tab_anim, text="Показать анимацию")
        self.anim_display.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.create_entry(self.tab_anim, "anim_int", "Интервал (мс):", 50, 0, col=1)

    def init_dxf_tab(self):
        self.dxf_import = ctk.CTkCheckBox(self.tab_dxf, text="Экспортировать в DXF", command=self.toggle_dxf_name)
        self.dxf_import.grid(row=0, column=0, columnspan=2, padx=10, pady=20, sticky="w")

        self.create_entry(self.tab_dxf, "dxf_name", "Имя файла DXF:", "kulachok_1", 1)

        ctk.CTkLabel(self.tab_dxf, text="Тип геометрии:").grid(row=2, column=0, padx=10, pady=5, sticky="w")
        self.dxf_line = ctk.CTkOptionMenu(self.tab_dxf, values=["spline", "line"])
        self.dxf_line.grid(row=2, column=1, padx=10, pady=5, sticky="w")

    def get_current_data_dict(self) -> Dict[str, Any]:
        """
        Собирает данные из виджетов и формирует словарь,
        структурно совпадающий с CamSolveOptions (и вложенными моделями).
        """
        calc_type = self.geom_calc_type.get()

        # Данные для CamConfig
        cam_config_data = {
            "N_k": float(self.entries["k_Nk"].get()),
            "D": float(self.entries["k_D"].get()),
            "h": float(self.entries["k_h"].get()),
            "z": float(self.entries["k_z"].get()),
            "D_t": float(self.entries["k_Dt"].get()),
            "R_r": float(self.entries["k_Rr"].get()),
            "f_pod": float(self.entries["k_fpod"].get()),
            "f_v": float(self.entries["k_fv"].get()),
            "f_op": float(self.entries["k_fop"].get()),
            "f_z": float(self.entries["k_fz"].get()),
        }

        # Данные для CalculatorConfig
        calculator_config_data = {}
        if calc_type == 'polidain':
            calculator_config_data = {
                "m": int(self.entries["pd_m"].get()),
                "d": int(self.entries["pd_d"].get()),
                "k_1": int(self.entries["pd_k1"].get()),
                "k_2": int(self.entries["pd_k2"].get()),
                "k_3": int(self.entries["pd_k3"].get()),
                "k_4": int(self.entries["pd_k4"].get()),
            }
        elif calc_type == 'polinmail':
            bc_str = self.entries["pm_bc"].get()
            bc_list = [float(x.strip()) for x in bc_str.split(",") if x.strip()]
            calculator_config_data = {
                "config_1": {
                    "m": int(self.entries["pm_m"].get()),
                    "d": int(self.entries["pm_d"].get()),
                    "boundary_conditions": bc_list
                }
            }

        # Собираем все опции
        data = {
            "cam_geometry_options": {
                "calculator_type": calc_type,
                "kulachok_type": self.geom_kulachok_type.get(),
                "N": int(self.entries["geom_n"].get()),
                "initial_angle": float(self.entries["geom_angle"].get()),
                "calculate_optimal_initial_angle": bool(self.geom_opt_angle.get()),
                "cam_config": cam_config_data,
                "calculator_config": calculator_config_data
            },
            "plotter_options": {
                "graphs_tolkatel_flag": bool(self.plot_tolkatel.get()),
                "graphs_kulachok_flag": bool(self.plot_kulachok.get()),
                "graphs_argument_type": self.plot_arg_type.get(),
                "graphs_profile_flag": bool(self.plot_profile.get()),
                "profile_and_graphs_together_flag": bool(self.plot_together.get()),
            },
            "rotate_animation_options": {
                "display_animation_flag": bool(self.anim_display.get()),
                "animation_intarval": int(self.entries["anim_int"].get()),
            },
            "dxf_creator_options": {
                "import_dxf_flag": bool(self.dxf_import.get()),
                "dxf_profile_name": self.entries["dxf_name"].get(),
                "dxf_line_type": self.dxf_line.get(),
            }
        }
        return data

    def generate_config(self):
        try:
            data = self.get_current_data_dict()

            # --- Создание Pydantic моделей ---

            # 1. CamGeometryOptions
            geom_data = data["cam_geometry_options"]

            # Создаем вложенные модели для CamGeometryOptions
            cam_conf = KulachokConfig(**geom_data["cam_config"])

            calc_conf = None
            if geom_data["calculator_type"] == 'polidain':
                calc_conf = PolidainConfig(**geom_data["calculator_config"])
            elif geom_data["calculator_type"] == 'polinmail':
                 # PolinmailConfig ожидает LocalPolinmailConfig в config_1
                 local_conf_data = geom_data["calculator_config"]["config_1"]
                 local_conf = LocalPolinmailConfig(**local_conf_data)
                 calc_conf = PolinmailConfig(config_1=local_conf)

            cam_geom_opts = CamGeometryOptions(
                calculator_type=geom_data["calculator_type"],
                kulachok_type=geom_data["kulachok_type"],
                N=geom_data["N"],
                initial_angle=geom_data["initial_angle"],
                calculate_optimal_initial_angle=geom_data["calculate_optimal_initial_angle"],
                cam_config=cam_conf,
                calculator_config=calc_conf
            )

            # 2. PlotterOptions
            plot_opts = PlotterOptions(**data["plotter_options"])

            # 3. RotateAnimationOptions
            anim_opts = RotateAnimationOptions(**data["rotate_animation_options"])

            # 4. DxfCreatorOptions
            dxf_opts = DxfCreatorOptions(**data["dxf_creator_options"])

            # 5. Итоговый объект
            final_config = CamSolveOptions(
                cam_geometry_options=cam_geom_opts,
                plotter_options=plot_opts,
                rotate_animation_options=anim_opts,
                dxf_creator_options=dxf_opts
            )

            # Запуск расчета
            calculate_cam_solve(final_config)

            messagebox.showinfo("Успех", "Расчет успешно выполнен!")

        except ValueError as e:
            import traceback
            traceback.print_exc()
            messagebox.showerror("Ошибка ввода", f"Проверьте правильность заполнения числовых полей.\nПодробности: {e}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            messagebox.showerror("Ошибка", str(e))

    def save_to_json(self):
        try:
            data = self.get_current_data_dict()
            filepath = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
            if filepath:
                with open(filepath, 'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=4, ensure_ascii=False)
                messagebox.showinfo("Успех", "Конфигурация сохранена!")
        except Exception as e:
            messagebox.showerror("Ошибка сохранения", str(e))

    def load_from_json(self):
        filepath = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if not filepath:
            return

        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)

            # --- CamGeometryOptions ---
            geom = data.get("cam_geometry_options", {})

            # Общие настройки
            calc_type = geom.get("calculator_type", "polidain")
            self.geom_calc_type.set(calc_type)
            self.update_method_visibility(calc_type)

            kulachok_type = geom.get("kulachok_type", "thin")
            self.geom_kulachok_type.set(kulachok_type)
            self.update_follower_params(kulachok_type)

            if "N" in geom: self.entries["geom_n"].delete(0, 'end'); self.entries["geom_n"].insert(0, str(geom["N"]))
            if "initial_angle" in geom: self.entries["geom_angle"].delete(0, 'end'); self.entries["geom_angle"].insert(0, str(geom["initial_angle"]))

            if "calculate_optimal_initial_angle" in geom:
                if geom["calculate_optimal_initial_angle"]:
                    self.geom_opt_angle.select()
                else:
                    self.geom_opt_angle.deselect()
                self.toggle_initial_angle()

            # Cam Config
            kul_conf = geom.get("cam_config", {})
            mapping_cam = {
                "N_k": "k_Nk", "D": "k_D", "h": "k_h", "z": "k_z", "D_t": "k_Dt", "R_r": "k_Rr",
                "f_pod": "k_fpod", "f_v": "k_fv", "f_op": "k_fop", "f_z": "k_fz"
            }
            for key, entry_key in mapping_cam.items():
                if key in kul_conf:
                    self.entries[entry_key].delete(0, 'end')
                    self.entries[entry_key].insert(0, str(kul_conf[key]))

            # Calculator Config
            calc_conf = geom.get("calculator_config", {})
            if calc_type == 'polidain':
                mapping_poly = {
                    "m": "pd_m", "d": "pd_d", "k_1": "pd_k1", "k_2": "pd_k2", "k_3": "pd_k3", "k_4": "pd_k4"
                }
                for key, entry_key in mapping_poly.items():
                    if key in calc_conf:
                        self.entries[entry_key].delete(0, 'end')
                        self.entries[entry_key].insert(0, str(calc_conf[key]))
            elif calc_type == 'polinmail':
                conf_1 = calc_conf.get("config_1", {})
                if "m" in conf_1: self.entries["pm_m"].delete(0, 'end'); self.entries["pm_m"].insert(0, str(conf_1["m"]))
                if "d" in conf_1: self.entries["pm_d"].delete(0, 'end'); self.entries["pm_d"].insert(0, str(conf_1["d"]))
                if "boundary_conditions" in conf_1:
                    bc_list = conf_1["boundary_conditions"]
                    bc_str = ", ".join(map(str, bc_list))
                    self.entries["pm_bc"].delete(0, 'end')
                    self.entries["pm_bc"].insert(0, bc_str)

            # --- PlotterOptions ---
            plot = data.get("plotter_options", {})
            if plot.get("graphs_tolkatel_flag"): self.plot_tolkatel.select()
            else: self.plot_tolkatel.deselect()

            if plot.get("graphs_kulachok_flag"): self.plot_kulachok.select()
            else: self.plot_kulachok.deselect()

            if plot.get("graphs_profile_flag"): self.plot_profile.select()
            else: self.plot_profile.deselect()

            if plot.get("profile_and_graphs_together_flag"): self.plot_together.select()
            else: self.plot_together.deselect()

            if "graphs_argument_type" in plot: self.plot_arg_type.set(plot["graphs_argument_type"])

            # --- RotateAnimationOptions ---
            anim = data.get("rotate_animation_options", {})
            if anim.get("display_animation_flag"): self.anim_display.select()
            else: self.anim_display.deselect()

            if "animation_intarval" in anim:
                self.entries["anim_int"].delete(0, 'end')
                self.entries["anim_int"].insert(0, str(anim["animation_intarval"]))

            # --- DxfCreatorOptions ---
            dxf = data.get("dxf_creator_options", {})
            if dxf.get("import_dxf_flag"):
                self.dxf_import.select()
            else:
                self.dxf_import.deselect()
            self.toggle_dxf_name()

            if "dxf_profile_name" in dxf:
                self.entries["dxf_name"].delete(0, 'end')
                self.entries["dxf_name"].insert(0, str(dxf["dxf_profile_name"]))

            if "dxf_line_type" in dxf:
                self.dxf_line.set(dxf["dxf_line_type"])

            messagebox.showinfo("Успех", "Данные успешно загружены!")
        except Exception as e:
            messagebox.showerror("Ошибка загрузки", str(e))


if __name__ == "__main__":
    app = CamConfiguratorApp()
    app.mainloop()
