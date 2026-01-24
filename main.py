import pysam
from tqdm import tqdm
import numpy as np
import math
from PIL import Image
import matplotlib.pyplot as plt

SLIDE_WIDTH = 3350
SLIDE_HEIGHT = 3350


def get_color(x: int, y: int):
    r = int(255 * x / SLIDE_WIDTH)
    g = int(255 * y / SLIDE_HEIGHT)
    b = int(255 * (y + x) / (SLIDE_WIDTH + SLIDE_HEIGHT))

    return r, g, b


def get_color_hue_based(x, y, width=100, height=100, cell_size=100):
    """
    Génère une couleur HSV unique pour chaque case d'une grille.
    Avec effet de vignette : clair au centre, sombre sur les bords.
    Maximise la diversité des couleurs entre les cases.

    Args:
        x: Coordonnée x
        y: Coordonnée y
        width: Largeur totale du carré (défaut: 100)
        height: Hauteur totale du carré (défaut: 100)
        cell_size: Taille d'une case de la grille (défaut: 100)

    Returns:
        Tuple (r, g, b) avec des valeurs entre 0 et 255
    """
    # Déterminer dans quelle case on se trouve
    cell_x = int(x / cell_size)
    cell_y = int(y / cell_size)

    # Générer une teinte unique pour cette case avec meilleure distribution
    # Utiliser plusieurs nombres premiers pour maximiser la variation
    seed1 = (cell_x * 73856093 + cell_y * 19349663) % 999983
    seed2 = (cell_x * 50331653 + cell_y * 83492791) % 999979
    seed3 = (cell_x * 12582917 + cell_y * 25165843) % 999961

    # Teinte : utiliser golden ratio pour une meilleure distribution des couleurs
    golden_ratio = 0.618033988749895
    hue = (seed1 * golden_ratio) % 1.0  # Entre 0 et 1, bien distribuée

    # Saturation variable pour plus de diversité (entre 0.6 et 1.0)
    saturation = 0.1 + (seed2 % 400) / 1000.0

    # Valeur de base pour la case (entre 0.7 et 0.95)
    base_value = 0.9 + (seed3 % 250) / 1000.0

    # Calculer la distance du pixel au centre du carré pour l'effet vignette
    center_x = width / 2
    center_y = height / 2
    dx = (x - center_x) / center_x
    dy = (y - center_y) / center_y
    distance_from_center = math.sqrt(dx * dx + dy * dy)

    # Ajuster la valeur (luminosité) en fonction de la distance
    vignette_factor = max(0.1, 1.0 - distance_from_center * 0.7)
    value = base_value * vignette_factor

    # Convertir HSV vers RGB
    def hsv_to_rgb(h, s, v):
        if s == 0:
            return (int(v * 255), int(v * 255), int(v * 255))

        i = int(h * 6.0)
        f = (h * 6.0) - i
        p = v * (1.0 - s)
        q = v * (1.0 - s * f)
        t = v * (1.0 - s * (1.0 - f))
        i = i % 6

        if i == 0:
            r, g, b = v, t, p
        elif i == 1:
            r, g, b = q, v, p
        elif i == 2:
            r, g, b = p, v, t
        elif i == 3:
            r, g, b = p, q, v
        elif i == 4:
            r, g, b = t, p, v
        else:
            r, g, b = v, p, q

        return (int(r * 255), int(g * 255), int(b * 255))

    return hsv_to_rgb(hue, saturation, value)


def get_color_checkerboard(x, y, width=100, height=100, cell_size=100):
    """
    Version avec un motif en damier où les cases alternent entre des couleurs
    très contrastées pour maximiser la distinction.
    """
    cell_x = int(x / cell_size)
    cell_y = int(y / cell_size)

    # Utiliser un pattern de damier pour alterner entre palettes
    is_even = (cell_x + cell_y) % 2 == 0

    # Seed pour la couleur
    seed = cell_x * 73856093 + cell_y * 19349663

    if is_even:
        # Couleurs chaudes
        hue = ((seed % 60) + 0) / 360.0  # Rouge à jaune
    else:
        # Couleurs froides
        hue = ((seed % 60) + 180) / 360.0  # Cyan à bleu

    saturation = 0.8 + (seed % 20) / 100.0
    value = 0.85 + (seed % 15) / 100.0

    # Convertir HSV vers RGB (même fonction que ci-dessus)
    def hsv_to_rgb(h, s, v):
        if s == 0:
            return (int(v * 255), int(v * 255), int(v * 255))

        i = int(h * 6.0)
        f = (h * 6.0) - i
        p = v * (1.0 - s)
        q = v * (1.0 - s * f)
        t = v * (1.0 - s * (1.0 - f))
        i = i % 6

        if i == 0:
            r, g, b = v, t, p
        elif i == 1:
            r, g, b = q, v, p
        elif i == 2:
            r, g, b = p, v, t
        elif i == 3:
            r, g, b = p, q, v
        elif i == 4:
            r, g, b = t, p, v
        else:
            r, g, b = v, p, q

        return (int(r * 255), int(g * 255), int(b * 255))

    return hsv_to_rgb(hue, saturation, value)


def get_color_gradient(x, y, taille=3350):
    """
    Calcule une couleur RGB en fonction de la position dans un carré.

    Args:
        x: Coordonnée x du pixel (0 à taille-1)
        y: Coordonnée y du pixel (0 à taille-1)
        taille: Taille du carré (défaut: 3350)

    Returns:
        tuple: (r, g, b) avec des valeurs entre 0 et 255
    """
    # Centrer les coordonnées (-0.5 à 0.5)
    cx = (x - taille / 2) / (taille / 2)
    cy = (y - taille / 2) / (taille / 2)

    # Distance au centre (0 au centre, ~0.707 aux coins)
    distance = np.sqrt(cx**2 + cy**2)

    # Normaliser la distance (0 au centre, 1 aux coins)
    distance_norm = min(distance / np.sqrt(2), 1.0)

    # Assigner des couleurs aux 4 quadrants
    # Coin haut-gauche: Rouge
    # Coin haut-droit: Vert
    # Coin bas-droit: Bleu
    # Coin bas-gauche: Magenta

    if cx <= 0 and cy <= 0:  # Haut-gauche
        base_color = np.array([255, 0, 0])  # Rouge
    elif cx > 0 and cy <= 0:  # Haut-droit
        base_color = np.array([0, 255, 0])  # Vert
    elif cx > 0 and cy > 0:  # Bas-droit
        base_color = np.array([0, 0, 255])  # Bleu
    else:  # Bas-gauche
        base_color = np.array([255, 0, 255])  # Magenta

    # Interpoler entre blanc (centre) et la couleur de base (coin)
    blanc = np.array([255, 255, 255])
    couleur = blanc * (1 - distance_norm) + base_color * distance_norm

    # Convertir en entiers
    r, g, b = couleur.astype(int)

    return (int(r), int(g), int(b))


def create_spatial_color(input_file: str, output_file: str):

    chrom = "chr22"
    start = 0
    end = 50000000

    with (
        pysam.AlignmentFile(input_file, "rb") as in_bam,
        pysam.AlignmentFile(output_file, "wb", header=in_bam.header) as out_bam,
    ):
        for read in tqdm(in_bam.fetch(chrom, start, end)):
            if not read.has_tag("CB"):
                continue

            barcode = read.get_tag("CB")

            _, _, x, y = barcode.split("_")
            y = y.replace("-1", "")
            x = int(x)
            y = int(y)

            r, g, b = get_color_hue_based(x, y, SLIDE_WIDTH, SLIDE_HEIGHT, SLIDE_WIDTH / 2)

            read.set_tag("YC", f"{r},{g},{b}")

            out_bam.write(read)


def get_color_a(x, y, width=100, height=100):
    """
    Calcule une couleur RGB basée sur la position (x, y) dans un carré.
    Crée plusieurs régions de couleurs différentes avec des transitions douces.
    Les points proches ont des couleurs similaires, les points éloignés ont des couleurs très différentes.

    Args:
        x: Coordonnée x (0 à width)
        y: Coordonnée y (0 à height)
        width: Largeur du carré (défaut: 100)
        height: Hauteur du carré (défaut: 100)

    Returns:
        Tuple (r, g, b) avec des valeurs entre 0 et 255
    """
    # Normaliser les coordonnées entre 0 et 1
    nx = x / width
    ny = y / height

    # Utiliser plusieurs fréquences pour créer des régions distinctes
    # Rouge: combinaison de 2 ondes sinusoïdales
    r = int(
        ((math.sin(nx * math.pi * 4) + 1) * 0.4 + (math.sin(ny * math.pi * 3) + 1) * 0.35) * 255
    )

    # Vert: combinaison différente
    g = int(
        ((math.sin(ny * math.pi * 5) + 1) * 0.4 + (math.sin((nx + ny) * math.pi * 2) + 1) * 0.35)
        * 255
    )

    # Bleu: motif en spirale
    distance = math.sqrt((nx - 0.5) ** 2 + (ny - 0.5) ** 2)
    angle = math.atan2(ny - 0.5, nx - 0.5)
    b = int(((math.sin(distance * math.pi * 8) + 1) * 0.4 + (math.sin(angle * 3) + 1) * 0.35) * 255)

    # S'assurer que les valeurs sont dans la plage valide
    r = max(0, min(255, r))
    g = max(0, min(255, g))
    b = max(0, min(255, b))

    return (r, g, b)


def get_color_perlin_like(x, y, width=100, height=100):
    """
    Version avec un effet type Perlin noise pour des régions organiques.
    """
    nx = x / width
    ny = y / height

    # Plusieurs octaves de bruit sinusoïdal
    def multi_sine(x, y, freq):
        return (
            math.sin(x * freq) * 0.5
            + math.sin(y * freq) * 0.5
            + math.sin((x + y) * freq * 0.7) * 0.3
        )

    # Rouge: basses fréquences
    r = int(((multi_sine(nx, ny, math.pi * 3) + 1.65) / 3.3) * 255)

    # Vert: fréquences moyennes
    g = int(((multi_sine(nx, ny, math.pi * 5) + 1.65) / 3.3) * 255)

    # Bleu: hautes fréquences
    b = int(((multi_sine(nx, ny, math.pi * 7) + 1.65) / 3.3) * 255)

    r = max(0, min(255, r))
    g = max(0, min(255, g))
    b = max(0, min(255, b))

    return (r, g, b)


def get_color_voronoi_like(x, y, width=100, height=100, num_regions=12):
    """
    Version inspirée des diagrammes de Voronoi pour créer des régions distinctes
    avec des transitions douces entre elles.
    """
    nx = x / width
    ny = y / height

    # Créer des "centres" de régions de manière déterministe
    regions = []
    for i in range(num_regions):
        # Utiliser une grille pseudo-aléatoire basée sur l'index
        cx = math.sin(i * 12.9898) * 0.5 + 0.5
        cy = math.sin(i * 78.233) * 0.5 + 0.5
        regions.append((cx, cy))

    # Calculer les distances à tous les centres avec un facteur de lissage
    weights = []
    total_weight = 0
    smoothing = 8.0  # Plus élevé = transitions plus douces

    for i, (cx, cy) in enumerate(regions):
        dist = math.sqrt((nx - cx) ** 2 + (ny - cy) ** 2)
        weight = math.exp(-dist * smoothing)
        weights.append(weight)
        total_weight += weight

    # Calculer la couleur comme moyenne pondérée
    r, g, b = 0, 0, 0
    for i, weight in enumerate(weights):
        normalized_weight = weight / total_weight
        # Couleur unique pour chaque région basée sur son index
        region_r = int((math.sin(i * 2.5) + 1) * 127.5)
        region_g = int((math.sin(i * 3.7 + 2) + 1) * 127.5)
        region_b = int((math.sin(i * 5.1 + 4) + 1) * 127.5)

        r += region_r * normalized_weight
        g += region_g * normalized_weight
        b += region_b * normalized_weight

    return (int(r), int(g), int(b))


def create_square():

    w = SLIDE_WIDTH
    h = SLIDE_HEIGHT

    img = Image.new("RGB", (w, h))

    for x in range(w):
        for y in range(h):
            img.putpixel((x, y), get_color_hue_based(x, y, w, h, w / 4))

    plt.imsave("square_grad.png", img)


if __name__ == "__main__":

    create_spatial_color(
        "../../dataset/pancreas/my_output/outs/possorted_genome_bam.bam", "test.bam"
    )

    # create_square()
