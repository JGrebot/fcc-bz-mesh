from email import parser
import numpy as np
from argparse import ArgumentParser


def get_path(path_points: str) -> list:
    """
    Get the path from the string of points
    """
    BZ_points = {
        "G": np.array([0, 0, 0]),
        "X":  np.array([0, 1/2, 0]),
        "L":  np.array([1/4, 1/4, 1/4]),
        "W":  np.array([1/4, 1/2, 0]),
        "U":  np.array([1/8, 1/2, 1/8]),
        "K":  np.array([3/8, 3/8, 0]),
    }
    path = []
    for point in list(path_points):
        path.append(BZ_points[point])
    return path


def build_path(path: list, n_points: int) -> np.array:
    """
    Build the path from the list of points
    """
    list_distance = []
    for ii in range(len(path) - 1):
        list_distance.append(np.linalg.norm(path[ii + 1] - path[ii]))
    list_distance = np.array(list_distance)
    total_distance = np.sum(list_distance)

    list_points = []
    for index_point, point in enumerate(path[:-1]):
        nb_points = int(n_points * list_distance[index_point] / total_distance)
        list_points += list(np.linspace(point,
                            path[index_point + 1], nb_points))
    return np.array(list_points)


def main():
    parser = ArgumentParser()
    parser.add_argument("-p", "--path", dest="path", type=str, required=True,
                        help="Path of the BZ. List of possible points:  (G, X, L, W, U, K)")
    parser.add_argument("-nb", "--nb", dest="nb", type=int, required=True)
    parser.add_argument("-o", "--output", dest="output",
                        type=str, required=False, default="path.dat")
    args = parser.parse_args()

    PATH_POINTS = get_path(args.path)
    PATH = build_path(PATH_POINTS, args.nb)
    list_x = PATH[:, 0]
    list_y = PATH[:, 1]
    list_z = PATH[:, 2]

    output_file = args.output
    print(f"Saving path in {output_file}")
    np.savetxt(output_file, np.transpose([list_x, list_y, list_z]), fmt="%.6f")


if __name__ == '__main__':
    main()
