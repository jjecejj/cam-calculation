from pydantic import BaseModel, Field
from typing import Literal

from core.cam_geometry import Kulachok
from exporters.dxf_creator.logic import build_profile


class DxfCreatorOptions(BaseModel):
    import_dxf_flag: bool = Field(default=False, description="Экспортировать профиль в формат DXF")
    dxf_profile_name: str = Field(default="kulachok_1", description="Имя файла DXF")
    dxf_line_type: Literal["spline", "line"] = Field(default="spline", description="Тип геометрии в DXF")

def resolve_dxf_creator_options(config: DxfCreatorOptions, kulachok: Kulachok):
    if config.import_dxf_flag:
        build_profile(
            kulachok.profile_data,
            profile_name=config.dxf_profile_name,
            line_type=config.dxf_line_type
        )