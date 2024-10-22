"""Data structures to manage information."""
from typing import Optional
from dataclasses import dataclass, asdict, InitVar
from pathlib import Path

@dataclass
class Environ():
    """Structures the go-continuum environment."""
    basedir: Path = Path('./')
    """Base directory."""
    uvdata: Optional[Path] = None
    """Directory where the generated MSs will be."""
    continuum: Optional[Path] = None
    """Directory for continuum images."""
    plots: Optional[Path] = None
    """Directory for plots."""
    check_env: InitVar[bool] = False
    """Check that directories exist."""

    def __post_init__(self, check_env):
        if self.uvdata is None:
            self.uvdata = self.basedir / 'uvdata'
        if self.plots is None:
            self.plots = self.basedir / 'plots'
        if check_env:
            self.do_check_env()

    def __getitem__(self, key):
        dict_form = asdict(self)
        return dict_form[key]

    def do_check_env(self):
        """Check that all directories exist."""
        for path in asdict(self).values():
            path.mkdir(exist_ok=True)

    def mkdir(self, name: str):
        """Make directory if needed."""
        self[name].mkdir(exist_ok=True)

@dataclass
class GoCoEnviron(Environ):
    """Structures the go-continuum environment."""
    basedir: Path = Path('./')
    """Base directory."""
    uvdata: Optional[Path] = None
    """Directory where the generated MSs will be."""
    dirty: Optional[Path] = None
    """Directory for dirty images."""
    continuum_control: Optional[Path] = None
    """Directory for continuum control images."""
    continuum: Optional[Path] = None
    """Directory for continuum images."""
    cubes: Optional[Path] = None
    """Directory for data cubes."""
    plots: Optional[Path] = None
    """Directory for plots."""
    check_env: InitVar[bool] = False
    """Check that directories exist."""

    def __post_init__(self, check_env):
        if self.dirty is None:
            self.dirty = self.basedir / 'dirty'
        if self.continuum_control is None:
            self.continuum_control = self.basedir / 'continuum_control'
        if self.continuum is None:
            self.continuum_auto = self.basedir / 'continuum_auto'
        if self.cubes is None:
            self.cubes = self.basedir / 'cubes'
        super().__post_init__(check_env)

