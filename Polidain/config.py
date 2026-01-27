from pathlib import Path
from dataclasses import dataclass, field
from matplotlib import rcParams
from core.schemas import PolidainConfig
import numpy as np

@dataclass
class PlotSettings:
    font_family: str
    font_size: int
    mathtext_fontset: str
    mathtext_rm: str
    mathtext_it: str
    rc: dict = field(default_factory=dict)
    font_serif: list = field(default_factory=list)

plotSettings_default = PlotSettings(
    rc = {"font.family": "serif", "mathtext.fontset": "stix"},
    font_serif = ["Times New Roman"] + rcParams["font.serif"],
    font_family = "Times New Roman",
    font_size = 12,
    mathtext_fontset = 'custom',
    mathtext_rm = "Times New Roman",
    mathtext_it = "Times New Roman:italic"
)

polidainConfig_default = PolidainConfig(
            N_k = 1000,
            D = 30.0 * 1e-3,
            D_t = 32.0 * 1e-3,
            h = 12.0 * 1e-3,
            z = 0.25 * 1e-3,
            f_pod = 80.0 / 180 * np.pi,
            f_v = 5.0 / 180 * np.pi,
            f_op = 75.0 / 180 * np.pi,
            f_z = 25 / 180 * np.pi,
            m = 3,
            d = 12,
            k_1 = 20,
            k_2 = 20,
            k_3 = 20,
            k_4 = 20
)

@dataclass
class AppConfig:
    cam: PolidainConfig = field(default_factory=lambda: polidainConfig_default)
    plot: PlotSettings = field(default_factory=lambda: plotSettings_default)

    # Автоматическое создание путей
    base_dir: Path = Path(__file__).parent
    output_dir: Path = base_dir / "data" / "output"

    def __post_init__(self):
        # Автоматически создать папку output, если её нет
        self.output_dir.mkdir(parents=True, exist_ok=True)