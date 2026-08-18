#!/bin/bash

chip_preset_names_csv() {
    printf '%s\n' "50-50, 50-20, 100-20"
}

apply_chip_preset() {
    case "$1" in
        50-50)
            x_spots_number=50
            y_spots_number=50
            length_spot=50
            interval=50
            banksy_lambda=${banksy_lambda:-0.6}
            banksy_resolution=${banksy_resolution:-1.5}
            banksy_spatial_neighbors=${banksy_spatial_neighbors:-15}
            banksy_cluster_neighbors=${banksy_cluster_neighbors:-25}
            banksy_subcluster_min_parent_spots=${banksy_subcluster_min_parent_spots:-40}
            banksy_subcluster_min_spots=${banksy_subcluster_min_spots:-10}
            banksy_subcluster_spatial_neighbors=${banksy_subcluster_spatial_neighbors:-8}
            banksy_subcluster_cluster_neighbors=${banksy_subcluster_cluster_neighbors:-12}
            banksy_subcluster_min_differential_cell_types=${banksy_subcluster_min_differential_cell_types:-4}
            ;;
        50-20)
            x_spots_number=50
            y_spots_number=50
            length_spot=20
            interval=20
            banksy_lambda=${banksy_lambda:-0.7}
            banksy_resolution=${banksy_resolution:-1.2}
            banksy_spatial_neighbors=${banksy_spatial_neighbors:-20}
            banksy_cluster_neighbors=${banksy_cluster_neighbors:-30}
            banksy_subcluster_min_parent_spots=${banksy_subcluster_min_parent_spots:-60}
            banksy_subcluster_min_spots=${banksy_subcluster_min_spots:-15}
            banksy_subcluster_spatial_neighbors=${banksy_subcluster_spatial_neighbors:-10}
            banksy_subcluster_cluster_neighbors=${banksy_subcluster_cluster_neighbors:-15}
            banksy_subcluster_min_differential_cell_types=${banksy_subcluster_min_differential_cell_types:-5}
            ;;
        100-20)
            x_spots_number=100
            y_spots_number=100
            length_spot=20
            interval=20
            banksy_lambda=${banksy_lambda:-0.8}
            banksy_resolution=${banksy_resolution:-1.0}
            banksy_spatial_neighbors=${banksy_spatial_neighbors:-30}
            banksy_cluster_neighbors=${banksy_cluster_neighbors:-50}
            banksy_subcluster_min_parent_spots=${banksy_subcluster_min_parent_spots:-100}
            banksy_subcluster_min_spots=${banksy_subcluster_min_spots:-25}
            banksy_subcluster_spatial_neighbors=${banksy_subcluster_spatial_neighbors:-15}
            banksy_subcluster_cluster_neighbors=${banksy_subcluster_cluster_neighbors:-25}
            banksy_subcluster_min_differential_cell_types=${banksy_subcluster_min_differential_cell_types:-6}
            ;;
        *) return 1 ;;
    esac

    banksy_pca_components=${banksy_pca_components:-20}
    banksy_max_m=${banksy_max_m:-1}
    banksy_neighbor_decay=${banksy_neighbor_decay:-scaled_gaussian}
    banksy_subcluster=${banksy_subcluster:-True}
    banksy_subcluster_max_depth=${banksy_subcluster_max_depth:-1}
    banksy_subcluster_resolution=${banksy_subcluster_resolution:-0.5}
    banksy_subcluster_max_dominant_fraction=${banksy_subcluster_max_dominant_fraction:-0.9}
    banksy_subcluster_lambda=${banksy_subcluster_lambda:-}
    banksy_subcluster_pca_components=${banksy_subcluster_pca_components:-}
}
