from pathlib import Path
from dataclasses import dataclass, field
from multiprocessing import freeze_support

from core.cam_geometry._config import default_kulachok_config
from core.profiling_methods.polidain._config import default_polidain_config
from vizualization.plotter._config import plotConfig_default, PlotConfig
from core.cam_geometry import Kulachok


class ValidationError(Exception):
    pass

@dataclass
class AppConfig:
    cam: Kulachok = field(default_factory=lambda: Kulachok(default_kulachok_config, default_polidain_config))
    plot: PlotConfig = field(default_factory=lambda: plotConfig_default)

    # Автоматическое создание путей
    base_dir: Path = Path(__file__).parent
    output_dir: Path = base_dir / "data" / "output"

    def __post_init__(self):
        # Автоматически создать папку output, если её нет
        self.output_dir.mkdir(parents=True, exist_ok=True)
        freeze_support()