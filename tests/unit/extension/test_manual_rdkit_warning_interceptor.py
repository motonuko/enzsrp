from rdkit.Chem import AllChem

from rxndata2.extension.rdkit_warning import RDKitWarningInterceptor
from tests.test_utils.test_default_path import TestDefaultPath


# should be run manually for log testing.
# If RDKitWarningInterceptor intercepts all logs you can not see the print statement.
# noinspection PyBroadException
def test_log2():
    interceptor = RDKitWarningInterceptor()
    try:
        with interceptor:
            AllChem.ReactionFromRxnFile(TestDefaultPath().test_data.joinpath('71500.rxn').as_posix())
    except Exception:
        pass
    print("can you see this text?")


# should be run manually
# def test_log():
#     interceptor = RDKitWarningInterceptor()
#     try:
#         with interceptor:
#             AllChem.ReactionFromRxnFile(DataPath2.test_data_for_testing.joinpath('71500.rxn').as_posix())
#     except Exception:
#         pass
#     raise ValueError("can you see this text?")


if __name__ == "__main__":
    test_log2()
    # test_log()
