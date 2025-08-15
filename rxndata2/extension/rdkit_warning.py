import sys
from io import StringIO
from typing import List, Type

from rdkit.rdBase import LogToPythonStderr


class RDKitWarningInterceptor:
    _instance = None
    _is_active = False  # Prevent the same class (even different instances) from being used in nested `with` statements.

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(RDKitWarningInterceptor, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        self._original_stderr = None
        self._sio = None
        self._supress = []
        LogToPythonStderr()

    def __enter__(self):
        if RDKitWarningInterceptor._is_active:
            raise RuntimeError("CustomClass is already in use within another 'with' block.")
        RDKitWarningInterceptor._is_active = True
        assert self._original_stderr is None and self._sio is None
        self._original_stderr = sys.stderr  # stderr is likely unique during runtime, but save it just in case it gets overwritten elsewhere.
        self._sio = sys.stderr = StringIO()  # StringIO likely avoids buffering.
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        assert self._original_stderr is not None and self._sio is not None, print(self._original_stderr, self._sio)
        try:
            self.check_error_immediately()
        except Exception as e:
            raise e
        finally:
            sys.stderr = self._original_stderr
            self._sio.close()
            self._refresh()
            RDKitWarningInterceptor._is_active = False

    def set_ignore_warning_exceptions(self, supress: List[Type[Exception]]):
        self._supress = supress

    def _pop_value(self) -> str:
        value = self._sio.getvalue()
        self._sio.truncate(0)
        self._sio.seek(0)
        return value

    def _refresh(self):
        self._original_stderr = None
        self._sio = None
        self._supress = []

    def check_error_immediately(self):
        value = self._pop_value()
        try:
            for line in value.split("\n"):
                if line == "":
                    continue
                if "Warning" not in line and "WARNING" not in line:
                    raise RuntimeError("unexpected")
                # NOTE: using 'if' instead of 'elif' because multiple warnings can be included in one line
                if "ambiguous stereochemistry" in value:
                    # Warning: ambiguous stereochemistry - linear bond arrangement  NOTE: more patterns?
                    raise AmbiguousStereoChemistryWarningException(value)
                if "molecule is tagged as 3D, but all Z coords are zero and 2D stereo markers have been found, marking the mol as 2D" in value:
                    raise Tagged3DBut2DMarkersFoundException(value)
                if "Proton(s) added/removed" in value:
                    raise ProtonsAddedOrRemovedWarningException(value)
                if "Metal was disconnected" in value:
                    raise MetalWasDisconnectedWarningException(value)
                if "Omitted undefined stereo" in value:
                    raise OmittedUndefinedStereoWarningException(value)
                if "Charges were rearranged" in value:
                    raise ChargesWereRearrangedWarningException(value)
                if "not removing hydrogen atom without neighbors" in value:
                    raise NotRemovingHydrogenAtomWithoutNeighbors(value)
                if "not removing hydrogen atom with dummy atom neighbors" in value:
                    raise NotRemovingHydrogenAtomWithDummyAtomNeighbors(value)

                raise ValueError(f"undefined warning has thrown! Add error handling. content: {value}")
        except Exception as e:
            if isinstance(e, tuple(self._supress)):
                print(f"\n LOG: The following warning is ignored \n {e}")
            else:
                raise e


class AmbiguousStereoChemistryWarningException(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class ProtonsAddedOrRemovedWarningException(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class MetalWasDisconnectedWarningException(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class OmittedUndefinedStereoWarningException(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class ChargesWereRearrangedWarningException(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class NotRemovingHydrogenAtomWithoutNeighbors(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class NotRemovingHydrogenAtomWithDummyAtomNeighbors(Exception):
    def __init__(self, message: str):
        super().__init__(message)


class Tagged3DBut2DMarkersFoundException(Exception):
    def __init__(self, message: str):
        super().__init__(message)
