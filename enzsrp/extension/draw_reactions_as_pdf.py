from typing import List, Optional

from PIL import ImageDraw
from rdkit.Chem import AllChem, Draw


def draw_reactions_as_pdf(smarts: List[str], labels: Optional[List[str]] = None):
    assert labels is None or len(smarts) == len(labels)
    images = []
    for i, smart in enumerate(smarts):
        reaction = AllChem.ReactionFromSmarts(smart)
        img = Draw.ReactionToImage(reaction)
        if labels is not None:
            draw = ImageDraw.Draw(img)
            draw.text((10, 10), labels[i], (0, 0, 0))
        images.append(img)

    images[0].save("reactions.pdf", save_all=True, append_images=images[1:])

#
# if __name__ == '__main__':
#     x = [
#         "C=CC(C)(C)c1[nH]c2ccccc2c1/C=C1/NC(=O)[C@H](C)NC1=O.CC(C)=CCOP(=O)([O-])OP(=O)([O-])[O-]>>C=CC(C)(C)c1[nH]c2c(CC=C(C)C)cccc2c1/C=C1/NC(=O)[C@H](C)NC1=O.O=P([O-])([O-])OP(=O)([O-])O",
#         "C=CC(C)(C)c1[nH]c2ccccc2c1/C=C1/NC(=O)[C@H](C)NC1=O.CC(C)=CCOP(=O)([O-])OP(=O)([O-])[O-]>>C=CC(C)(C)c1[nH]c2ccc(CC=C(C)C)cc2c1/C=C1/NC(=O)[C@H](C)NC1=O.O=P([O-])([O-])OP(=O)([O-])O"]
#     draw_reactions_as_pdf(x, ["1", "2"])
