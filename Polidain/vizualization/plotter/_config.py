from dataclasses import dataclass, field
from matplotlib import rcParams

@dataclass
class PlotConfig:
    font_family: str
    font_size: int
    mathtext_fontset: str
    mathtext_rm: str
    mathtext_it: str
    rc: dict = field(default_factory=dict)
    font_serif: list = field(default_factory=list)

plotConfig_default = PlotConfig(
    rc = {"font.family": "serif", "mathtext.fontset": "stix"},
    font_serif = ["Times New Roman"] + rcParams["font.serif"],
    font_family = "Times New Roman",
    font_size = 12,
    mathtext_fontset = 'custom',
    mathtext_rm = "Times New Roman",
    mathtext_it = "Times New Roman:italic"
)
