from pydantic import BaseModel, Field
from typing import Literal


class DxfCreatorOptions(BaseModel):
    import_dxf_flag: bool = Field(default=False, description="Экспортировать профиль в формат DXF")
    dxf_profile_name: str = Field(default="kulachok_1", description="Имя файла DXF")
    dxf_line_type: Literal["spline", "line"] = Field(default="spline", description="Тип геометрии в DXF")
