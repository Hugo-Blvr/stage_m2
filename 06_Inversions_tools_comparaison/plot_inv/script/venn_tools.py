import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Ellipse
from matplotlib.path import Path
import matplotlib.patches as patches

def venn_4(labels, region_data, title, ax):
    # Paramètres des ellipses - ajustés pour mieux faire apparaître toutes les intersections
    centers = [
        (0.40, 0.50),  # Centre de l'ellipse A 
        (0.48, 0.56),  # Centre de l'ellipse B 
        (0.54, 0.56),  # Centre de l'ellipse C 
        (0.62, 0.50)   # Centre de l'ellipse D 
    ]
    
    # Largeur, hauteur et angle de rotation pour chaque ellipse
    widths = [0.55, 0.55, 0.55, 0.55]  # Largeur
    heights = [0.30, 0.30, 0.30, 0.30]  # Hauteur
    angles = [130, 130, 50, 50]  # Angles de rotation
    colors = ['blue', 'green', 'red', 'purple']
    #colors = ["#3498db", "#2ecc71", "#e74c3c", "#9b59b6"]

    label_positions = [
    (0.12, 0.45),  # Position du label A
    (0.25, 0.82),  # Position du label B  
    (0.70, 0.82),  # Position du label C
    (0.91, 0.45)   # Position du label D
    ]
    
    # Création des ellipses
    ellipses = []
    for center, width, height, angle, color, label, label_pos in zip(centers, widths, heights, angles, colors, labels, label_positions):
        ellipse = Ellipse(center, width, height, angle=angle, alpha=0.4, color=color, label=label)
        ax.add_patch(ellipse)
        ellipses.append(ellipse)
                
        # Ajout du label à la position personnalisée
        #plt.text(label_pos[0], label_pos[1], label, fontsize=10, fontweight='bold', ha='center', va='center')
    # Dictionnaire pour stocker les coordonnées des régions
    region_centers = {
        # Régions individuelles (A, B, C, D)
        (0,): (0.25, 0.60),  # A uniquement 
        (1,): (0.38, 0.75),  # B uniquement 
        (2,): (0.64, 0.75),  # C uniquement 
        (3,): (0.77, 0.60),  # D uniquement 
                
        # Intersections de 2 ensembles
        (0, 1): (0.35, 0.65),  # A ∩ B 
        (0, 2): (0.37, 0.40),  # A ∩ C 
        (0, 3): (0.51, 0.30),  # A ∩ D 
        (1, 2): (0.51, 0.67),  # B ∩ C 
        (1, 3): (0.65, 0.40),  # B ∩ D 
        (2, 3): (0.67, 0.65),  # C ∩ D 
                
        # Intersections de 3 ensembles
        (0, 1, 2): (0.44, 0.60),  # A ∩ B ∩ C 
        (0, 1, 3): (0.57, 0.35),  # A ∩ B ∩ D 
        (0, 2, 3): (0.45, 0.35),  # A ∩ C ∩ D 
        (1, 2, 3): (0.58, 0.60),  # B ∩ C ∩ D 
                
        # Intersection des 4 ensembles
        (0, 1, 2, 3): (0.51, 0.47)  # A ∩ B ∩ C ∩ D 
    }


    # Ajouter les valeurs à chaque région
    for region, value in region_data.items():
        region_text = f"{value}"
        
        # Obtenir les coordonnées de la région
        if region in region_centers:
            x, y = region_centers[region]
            # Ajouter le texte avec un fond pour une meilleure lisibilité
            ax.text(x, y, region_text, 
                    ha='center', va='center', 
                    fontsize=8, weight='bold',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))
        else:print(f"Avertissement: La région {region} n'a pas de coordonnées définies.")

    # Configuration des axes
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(title) 
    ax.set_aspect('equal')
    ax.axis('off')

    return ax


if __name__ == "__main__":
    data = {
    (0,): "A",        # A
    (1,): "B",        # B 
    (2,): "C",        # C 
    (3,): "D",        # D
    
    (0, 1): "A ∩ B",       # A ∩ B
    (0, 2): "A ∩ C",       # A ∩ C
    (0, 3): "A ∩ D",       # A ∩ D
    (1, 2): "B ∩ C",       # B ∩ C
    (1, 3): "B ∩ D",       # B ∩ D
    (2, 3): "C ∩ D",       # C ∩ D
    
    (0, 1, 2): "A ∩ B ∩ C",    # A ∩ B ∩ C
    (0, 1, 3): "A ∩ B ∩ D",     # A ∩ B ∩ D
    (0, 2, 3): "A ∩ C ∩ D",     # A ∩ C ∩ D
    (1, 2, 3): "B ∩ C ∩ D",     # B ∩ C ∩ D
    
    (0, 1, 2, 3): "A ∩ B ∩ C ∩ D"   # A ∩ B ∩ C ∩ D
    }

