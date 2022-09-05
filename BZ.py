#!/usr/bin/env python3

from argparse import ArgumentParser
import gmsh
import sys
import numpy as np
from numpy import linalg as LA

model = gmsh.model
factory = model.occ

pi = np.pi
scalingFactor = 2.0

DESCRIPTION_SCRIPT = """
This function builds a mesh of face cubic centered BZ. Though mesh_type you can
select either 1, 1/2, 1/8 or 1/48 of the BZ. mesh_name controls the name of
the written mesh. You can visualize the geometry and the mesh though -vm
and -gm options. Then you can add refinements box with the --refine_L and
--refine_delta options.
"""


def main_BZ_generating_mesh():
    """
    This function builds a mesh of the BZ of Si. Though mesh_type you can
    select either 1, 1/2, 1/8 or 1/48 of the BZ. mesh_name control the name of
    the written mesh. You can visualize the geometry and the mesh though -vm
    and -gm options. Then you can add refinements box with the --refine_L and
    --refine_delta option.

    The easier way to launch this script and understand its arguments is to
    type directly: python3 BZ.py --help

    Description of arguments can be found with: python3 BZ.py --help
    """
    ##############################################
    # Parsing command lines arguments
    ##############################################
    parser = ArgumentParser(
        description=DESCRIPTION_SCRIPT)
    parser.add_argument("-t",
                        "--mesh_type",
                        dest="mesh_type",
                        type=int, required=True,
                        help="""Type of mesh to generate, either 1, 2, 8 or 48.
                        For full BZ, half of the bz and so on.""")
    parser.add_argument("-n",
                        "--mesh_name",
                        dest="filename",
                        default="bz",
                        required=True,
                        help="""Name of mesh file to write. mesh_type is always
                        append to it.""")
    parser.add_argument("-lmin",
                        "--lengthMin",
                        dest="lengthMin",
                        type=float,
                        default=3e-2,
                        required=True,
                        help="""Minimal characteristic length of mesh elements (
                        [Ga, X] is of length 1/2).""")
    parser.add_argument("-lmax",
                        "--lengthMax",
                        dest="lengthMax",
                        type=float,
                        default=10e-2,
                        required=True,
                        help="""Maximal characteristic length of mesh elements (
                        [Ga, X] is of length 1/2).""")
    parser.add_argument("-vm",
                        dest="view_mesh",
                        action="store_true",
                        default=False,
                        help="""Use this option to open gmsh gui in order to
                        view the mesh.""")
    parser.add_argument("-gm",
                        dest="generate_mesh",
                        action="store_false",
                        default=True,
                        help="""Use this option to NOT generate the mesh. With
                        -vm also activated you'll see only the geometrical
                        model.""")
    parser.add_argument("-e", "--extension",
                        dest="extension",
                        default=".msh",
                        help="""Extension of the mesh file saved. The extension must be supported by GMSH.""")
    parser.add_argument("--refine_delta",
                        dest="refine_delta",
                        action='append',
                        nargs='+',
                        help="""Refinement parameters for delta valley. Three
                        parameters are required: h_fine, dL and dH. The
                        refinement box is centered in 0.875*[Ga, X]. dL is the
                        length along the [Ga, X] segment. dL must be between
                        0.0 and 0.125. dH is relative to the segment [X, W],
                        and must be inside ]0.0, \\frac{1}{\srqt{2}}[  ~=
                        ]0.0, 0.7[.""")
    parser.add_argument("--refine_L",
                        dest="refine_L",
                        action='append',
                        nargs='+',
                        help="""Refinement parameters for L valley. Three 
                        parameters are required: h_fine, dL and dH. The
                        refinement box is centered in 0.875*[Ga, L]. dL must be
                        between 0.0 and 0.125. dH is relative to the segment
                        [L, U] and must be inside ]0.0, 1.0].""")

    args = parser.parse_args()

    mesh_type = f"{args.mesh_type}"
    filename = f"{args.filename}_{mesh_type}"
    view_mesh = args.view_mesh
    generate_mesh = args.generate_mesh
    lengthMin = args.lengthMin
    lengthMax = args.lengthMax
    mesh_extension = args.extension

    if args.refine_L is not None:
        h_fine_L = float(args.refine_L[0][0])
        dL_L = float(args.refine_L[0][1])
        dH_L = float(args.refine_L[0][2])

    if args.refine_delta is not None:
        print(args.refine_delta)
        h_fine_delta = float(args.refine_delta[0][0])
        dL_delta = float(args.refine_delta[0][1])
        dH_delta = float(args.refine_delta[0][2])

    ##############################################
    # Parameters examples
    ##############################################
    # mesh_type = "1"   # Half of the BZ
    # mesh_type = "2"   # Half of the BZ
    # mesh_type = "8"   # 1/8 of the BZ
    # mesh_type = "48"  # 1/48 of the BZ

    # filename = "bz" + "_" + mesh_type

    # Do you wish to see the mesh ?
    # view_mesh = True

    # Do you wish to generate the mesh ?
    # generate_mesh = True

    # Ga -> X is 1/2
    # lengthMin = 3e-2
    # lengthMax = 10e-2

    # Local refinement parameters
    # h_fine_L = 0.01*lengthMin
    # dL_L = 0.125
    # dH_L = 0.10
    # h_fine_delta = 0.01*lengthMin
    # dL_delta = 0.125
    # dH_delta = 0.10

    ##############################################
    # Geometrical model building
    ##############################################
    initialize_gmsh_model(lengthMin, lengthMax, scalingFactor)
    BZ_points = build_BZ_points()
    dimTag = build_1_over_48_of_the_BZ(BZ_points)
    apply_symmetry_and_rotation(BZ_points, mesh_type, dimTag)

    ##############################################
    # Applying refinements
    ##############################################
    if args.refine_delta is not None:
        refine_mesh_delta_valley(
            dH_delta, dL_delta, mesh_type, BZ_points, h_fine_delta)
    if args.refine_L is not None:
        refine_mesh_L_valley(dH_L, dL_L, mesh_type, BZ_points, h_fine_L)

    ##############################################
    # Mesh Writing
    ##############################################
    writeMeshFile(filename, view=view_mesh,
                  generate=generate_mesh, extension=mesh_extension)


def build_BZ_points():
    """
    This function builds all the 3d points we need on the Brillouin zone in
    order to build 1/48 of it. So it return, stored in a dict, the following:
    Ga, L, U, X, W, K.

    The ScalingFactor (set in initialize_gmsh_model) is 4*pi/a . So all the
    coordinates are multiplied by this scalingFactor.

    output:
    BZ_points       dict: All keys are Ga, L, U, X, W, K. 
                          And each key is np array of 3 coordinates (Kx, Ky, Kz).
    """

    # BZ_points in term of a and pi.
    # BZ_points = {
    #         "Ga": [0, 0, 0],
    #         "X": [0, 2*pi/a, 0],
    #         "L": [pi/a, pi/a, pi/a],
    #         "W": [pi/a, 2*pi/a, 0],
    #         "U": [pi/(2*a), 2*pi/a, pi/(2*a)],
    #         "K": [3*pi/(2*a), 3*pi/(2*a), 0],
    # }
    BZ_points = {
        "Ga": np.array([0, 0, 0]),
        "X":  np.array([0, 1/2, 0]),
        "L":  np.array([1/4, 1/4, 1/4]),
        "W":  np.array([1/4, 1/2, 0]),
        "U":  np.array([1/8, 1/2, 1/8]),
        "K":  np.array([3/8, 3/8, 0]),
    }
    return BZ_points


def apply_symmetry_and_rotation(BZ_points, mesh_type, dimTag):
    """
    This function applies the rotation and symmetries to 1/48 of the BZ.
    Three type of operations are used:
        0) If 1/48 is chosen, do nothing.
        1) From 1/48 to 1/8 of the BZ. A symmetry along Ga-L-K plane and
           3 rotations along the Gamma-L axis
        2) From 1/8 to 1/2 of the BZ, 4 rotations on the z axis.
        3) From 1/2 to 1 of the BZ, symmetry along the plane of eq. z=0 .

    input:
    BZ_points       dict: All keys are Ga, L, U, X, W, K. 
                          And each key is np array of 3 coordinates (Kx, Ky, Kz).
    mesh_type       str:  Either 1, 2, or 8 or 48 for the three types of mesh you wish to generate.
    dimTag          tuple: (dimension, tag) of the 1/48 BZ volume previously created.

    output:
    (None) Information is stored in gmsh.model.
    """
    if mesh_type == "48":
        return None
    elif mesh_type == "8" or mesh_type == "2" or mesh_type == "1":
        Ga = BZ_points["Ga"]
        L = BZ_points["L"]

        # symmetry by the Ga - L - K plane (of equation x - y = 0)
        dimTag_2 = factory.copy(dimTag)
        factory.mirror(dimTag_2, 1, -1, 0, 0)
        dimTag_24 = factory.fuse(
            dimTag_2, dimTag, removeObject=True, removeTool=True)[0]
        # All dim-tag are of format: [(3, tag)]

        # Rotation of axis Ga -> L of angle 120 and 240
        dimTag_24_copy = factory.copy(dimTag_24)
        factory.rotate(
            dimTag_24_copy, Ga[0], Ga[1], Ga[2], L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], 2*pi/3)
        dimTag_24_copy2 = factory.copy(dimTag_24)
        factory.rotate(
            dimTag_24_copy2, Ga[0], Ga[1], Ga[2], L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], -2*pi/3)

        # Fusing created volumes
        dimTag_12 = factory.fuse(
            dimTag_24_copy, dimTag_24, removeObject=True, removeTool=True)[0]
        dimTag_8 = factory.fuse(
            dimTag_24_copy2, dimTag_12, removeObject=True, removeTool=True)[0]

        if mesh_type == "2" or mesh_type == "1":

            # Rotation of z-axis of angles 90, 180, 270 degrees.
            dimTag_8_copy1 = factory.copy(dimTag_8)
            factory.rotate(dimTag_8_copy1, 0, 0, 0, 0, 0, 1, pi/2)
            dimTag_8_copy2 = factory.copy(dimTag_8)
            factory.rotate(dimTag_8_copy2, 0, 0, 0, 0, 0, 1, 2*pi/2)
            dimTag_8_copy3 = factory.copy(dimTag_8)
            factory.rotate(dimTag_8_copy3, 0, 0, 0, 0, 0, 1, 3*pi/2)

            # Fusing created volumes
            dimTag_4_fuse1 = factory.fuse(
                dimTag_8, dimTag_8_copy1, removeObject=True, removeTool=True)[0]
            dimTag_4_fuse2 = factory.fuse(
                dimTag_4_fuse1, dimTag_8_copy2, removeObject=True, removeTool=True)[0]
            dimTag_2 = factory.fuse(
                dimTag_4_fuse2, dimTag_8_copy3, removeObject=True, removeTool=True)[0]

            if mesh_type == "1":
                # Symmetry along the plane of eq. z = 0
                dimTag_2_copy = factory.copy(dimTag_2)
                factory.mirror(dimTag_2_copy, 0, 0, 1, 0)

                # Fusing created volumes
                dimTag_1 = factory.fuse(
                    dimTag_2_copy, dimTag_2, removeObject=True, removeTool=True)[0]

    factory.synchronize()
    return None


def refine_mesh_L_valley(dH, dL, mesh_type, BZ_points, h_fine):
    """
    This functions add refinement boxes in all the L valley. 
    It works for all mesh_type automatically.

    The box along the Ga -> L axis is centered at 0.875*(Ga-L). Its length
    along the (Ga,L) axis is (L-Ga)*dL. dL cannot be higher than 0.125*||L-Ga||. 
    If this is the case we automatically set dL to its maximum.

    The box has size (K-L)*dH in both (L,K) and (L,U) axis.

    Said otherwise, the input dL and dH are normalized quantities.

    input:
    dH          double: box height (relative) along z and x axis.
    dL          double: box length (relative ) along y axis (X-Ga). Cannot be longer to (1-0.875)*||Ga-X||.
                        The length is automatically adapted.
    mesh_type   str: either 1, 2, or 8 or 48 for the three types of mesh you wish to generate.
    BZ_points   dict: All keys are Ga, L, U, X, W, K. 
                      And each key is np array of 3 coordinates (Kx, Ky, Kz).
    h_fine    double: Value of max tetrahedra size. Must be coherent with the parameter
                        lengthMin and lengthMax.

    output:     None (stored in gmsh.model)
    """
    Ga = BZ_points["Ga"]
    W = BZ_points["W"]
    U = BZ_points["U"]
    K = BZ_points["K"]
    L = BZ_points["L"]

    # Maximum value for dL (the refinement box is always centered in
    # 0.875(L-Ga).
    if dL > 0.125:
        dL = 0.125

    if dH > 1:
        dH = 1
    elif dH <= 0:
        print(f"{dH =}")
        print("Please read docstring and put a value for dH between 0 and 1.")
        sys.exit()
        printPanda()

    if mesh_type == "48" or mesh_type == "8" or mesh_type == "2" or mesh_type == "1":
        # Only one box around the L points. ORDER MATTER !
        all_points = []
        p0 = (L-Ga)*(0.875 - dL)
        p1 = np.copy(p0) + (U-L)*dH
        p2 = (L-Ga)*(0.875 + dL) + (U-L)*dH
        p3 = (L-Ga)*(0.875 + dL)
        p4 = np.copy(p3) + (K-L)*dH
        p5 = np.copy(p4) - (L-Ga)*2*dL
        all_points = [p0, p1, p2, p3, p4, p5]

        # Adding points. ORDER MATTER !
        tag_point = []
        for p in all_points:
            tag_point.append(factory.addPoint(p[0], p[1], p[2]))

        # Adding line. ORDER MATTER !
        tag_line = []

        tag_line.append(factory.addLine(tag_point[0], tag_point[1]))  # 0
        tag_line.append(factory.addLine(tag_point[1], tag_point[2]))  # 1
        tag_line.append(factory.addLine(tag_point[2], tag_point[3]))  # 2
        tag_line.append(factory.addLine(tag_point[3], tag_point[0]))  # 3

        tag_line.append(factory.addLine(tag_point[3], tag_point[4]))  # 4
        tag_line.append(factory.addLine(tag_point[4], tag_point[5]))  # 5
        tag_line.append(factory.addLine(tag_point[5], tag_point[0]))  # 6

        tag_line.append(factory.addLine(tag_point[1], tag_point[5]))  # 7
        tag_line.append(factory.addLine(tag_point[2], tag_point[4]))  # 8

        factory.synchronize()

        # Adding curve loop.
        tag_curveloop = []

        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[0], tag_line[1], tag_line[2], tag_line[3]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[2], tag_line[4], -tag_line[8]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[0], tag_line[7], -tag_line[6]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[1], tag_line[8], tag_line[5], -tag_line[7]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[3], -tag_line[6], -tag_line[5], -tag_line[4]]))

        # Adding Surface
        tag_surface = []
        for l in tag_curveloop:
            tag_surface.append(factory.addPlaneSurface([l]))

        # Adding Surface loop
        tag_surfaceloop = factory.addSurfaceLoop(tag_surface)

        # Adding Volume
        Tag_refinement = factory.addVolume([tag_surfaceloop])
        dimTag = [(3, Tag_refinement)]

        if mesh_type == "48":
            computeIntersectingVolumes()
            setSize(Tag_refinement, h_fine)

        if mesh_type == "8" or mesh_type == "2" or mesh_type == "1":
            Ga = BZ_points["Ga"]
            L = BZ_points["L"]

            # symmetry by the Ga - L - K plane (of equation x - y = 0)
            dimTag_2 = factory.copy(dimTag)
            factory.mirror(dimTag_2, 1, -1, 0, 0)
            dimTag_24 = factory.fuse(
                dimTag_2, dimTag, removeObject=True, removeTool=True)[0]

            # Rotation of axis Ga -> L of angle 120 and 240
            dimTag_3 = factory.copy(dimTag_24)
            dimTag_4 = factory.copy(dimTag_24)
            factory.rotate(dimTag_3, Ga[0], Ga[1], Ga[2],
                           L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], -2*pi/3)
            factory.rotate(dimTag_4, Ga[0], Ga[1], Ga[2],
                           L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], 2*pi/3)

            dimTag_24_tmp = factory.fuse(
                dimTag_24, dimTag_3, removeObject=True, removeTool=True)[0]
            dimTag_8 = factory.fuse(
                dimTag_24_tmp, dimTag_4, removeObject=True, removeTool=True)[0]

            if mesh_type == "8":
                computeIntersectingVolumes()
                setSize(dimTag_8[0][1], h_fine)

            if mesh_type == "2" or mesh_type == "1":
                # dimTag_8 is the only tag we have
                #
                # The rotation occurs along the z axis.
                #
                # We want to copy and rotate 3 time.
                # No fuse.

                # dimTag_8
                dimTag_8_1 = factory.copy(dimTag_8)
                dimTag_8_2 = factory.copy(dimTag_8)
                dimTag_8_3 = factory.copy(dimTag_8)
                factory.rotate(dimTag_8_1, 0, 0, 0, 0, 0, 1, pi/2)
                factory.rotate(dimTag_8_2, 0, 0, 0, 0, 0, 1, 2*pi/2)
                factory.rotate(dimTag_8_3, 0, 0, 0, 0, 0, 1, 3*pi/2)

                if mesh_type == "2":
                    computeIntersectingVolumes()
                    setSize(dimTag_8[0][1],   h_fine)
                    setSize(dimTag_8_1[0][1], h_fine)
                    setSize(dimTag_8_2[0][1], h_fine)
                    setSize(dimTag_8_3[0][1], h_fine)

                if mesh_type == "1":
                    # Mirror by the plane of equation z = 0
                    dimTag_2 = factory.copy(dimTag_8)
                    dimTag_2_1 = factory.copy(dimTag_8_1)
                    dimTag_2_2 = factory.copy(dimTag_8_2)
                    dimTag_2_3 = factory.copy(dimTag_8_3)

                    factory.mirror(dimTag_2,   0, 0, 1, 0)
                    factory.mirror(dimTag_2_1, 0, 0, 1, 0)
                    factory.mirror(dimTag_2_2, 0, 0, 1, 0)
                    factory.mirror(dimTag_2_3, 0, 0, 1, 0)

                    computeIntersectingVolumes()

                    setSize(dimTag_2[0][1], h_fine)
                    setSize(dimTag_2_1[0][1], h_fine)
                    setSize(dimTag_2_2[0][1], h_fine)
                    setSize(dimTag_2_3[0][1], h_fine)

                    setSize(dimTag_8[0][1], h_fine)
                    setSize(dimTag_8_1[0][1], h_fine)
                    setSize(dimTag_8_2[0][1], h_fine)
                    setSize(dimTag_8_3[0][1], h_fine)


def refine_mesh_delta_valley(dH, dL, mesh_type, BZ_points, h_fine):
    """
    This functions add refinement boxes in all the delta valley. 
    It works for all mesh_type automatically.

    The box along the Ga -> X axis is centered at 0.875 (Ga-X). Its length
    along the X- Ga axis is (X-Ga)*dL. dL cannot be higher than 0.125*||Ga-X||.
    If this is the case, we automatically set dL to its maximum.

    The box has size dH in both x-axis and z-axis.
    Rq: Ga -> X is colinear with y-axis.

    input:
    dH          double: box height along z and x axis.
    dL          double: box length (relative) along y axis (X-Ga). Cannot be longer to (1-0.875)*||Ga-X||.
                        The length is automatically adapted.
    mesh_type   str: either 1, 2, or 8 or 48 for the three types of mesh you wish to generate.
    BZ_points   dict: All keys are Ga, L, U, X, W, K. 
                      And each key is np array of 3 coordinates (Kx, Ky, Kz).
    h_fine    double: Value of max tetrahedra size. Must be coherent with the parameter
                        lengthMin and lengthMax.

    output:     None (stored in gmsh.model)
    """
    Ga = BZ_points["Ga"]
    X = BZ_points["X"]
    W = BZ_points["W"]
    U = BZ_points["U"]

    # Maximum value for dL (the refinement box is always centered in
    # 0.875(X-Ga).
    if dL > 0.125:
        dL = 0.125

    if dH > 1:
        dH = 1
    elif dH <= 0:
        print(f"{dH =}")
        print("Please read docstring and put a value for dH between 0 and 1.")
        sys.exit()
        printPanda()

    if mesh_type == "48" or mesh_type == "8" or mesh_type == "2" or mesh_type == "1":
        # Only one box around the X points.
        p0 = (X-Ga)*(0.875 - dL)
        p1 = np.copy(p0) + (W-X)*dH
        p2 = (X-Ga)*(0.875 + dL) + (W-X)*dH
        p3 = (X-Ga)*(0.875 + dL)
        p4 = (X-Ga)*(0.875 + dL) + (U-X)/LA.norm(U-X) * \
            LA.norm(W-X)*np.sqrt(dH**2 + dH**2)
        p5 = np.copy(p4) - (X-Ga)*(2*dL)
        all_points = [p0, p1, p2, p3, p4, p5]

        # Adding points. ORDER MATTER !
        tag_point = []
        for p in all_points:
            tag_point.append(factory.addPoint(p[0], p[1], p[2]))

        # Adding line. ORDER MATTER !
        tag_line = []

        tag_line.append(factory.addLine(tag_point[0], tag_point[1]))  # 0
        tag_line.append(factory.addLine(tag_point[1], tag_point[2]))  # 1
        tag_line.append(factory.addLine(tag_point[2], tag_point[3]))  # 2
        tag_line.append(factory.addLine(tag_point[3], tag_point[0]))  # 3

        tag_line.append(factory.addLine(tag_point[3], tag_point[4]))  # 4
        tag_line.append(factory.addLine(tag_point[4], tag_point[5]))  # 5
        tag_line.append(factory.addLine(tag_point[5], tag_point[0]))  # 6

        tag_line.append(factory.addLine(tag_point[1], tag_point[5]))  # 7
        tag_line.append(factory.addLine(tag_point[2], tag_point[4]))  # 8

        factory.synchronize()

        # Adding curve loop.
        tag_curveloop = []

        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[0], tag_line[1], tag_line[2], tag_line[3]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[2], tag_line[4], -tag_line[8]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[0], tag_line[7], -tag_line[6]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[1], tag_line[8], tag_line[5], -tag_line[7]]))
        tag_curveloop.append(factory.addCurveLoop(
            [tag_line[3], -tag_line[6], -tag_line[5], -tag_line[4]]))

        # Adding Surface
        tag_surface = []
        for l in tag_curveloop:
            tag_surface.append(factory.addPlaneSurface([l]))

        # Adding Surface loop
        tag_surfaceloop = factory.addSurfaceLoop(tag_surface)

        # Adding Volume
        Tag_refinement = factory.addVolume([tag_surfaceloop])
        dimTag = [(3, Tag_refinement)]

        if mesh_type == "48":
            computeIntersectingVolumes()
            setSize(Tag_refinement, h_fine)

        if mesh_type == "8" or mesh_type == "2" or mesh_type == "1":

            Ga = BZ_points["Ga"]
            L = BZ_points["L"]

            # symmetry by the Ga - L - K plane (of equation x - y = 0)
            dimTag_2 = factory.copy(dimTag)  # [(3, tag)]
            factory.mirror(dimTag_2, 1, -1, 0, 0)
            dimTag_24_1 = factory.fuse(
                dimTag_2, dimTag, removeObject=True, removeTool=True)[0]

            # Rotation of axis Ga -> L of angle 120 and 240
            dimTag_3 = factory.copy(dimTag)
            dimTag_5 = factory.copy(dimTag)
            factory.rotate(dimTag_3, Ga[0], Ga[1], Ga[2],
                           L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], -2*pi/3)
            factory.rotate(dimTag_5, Ga[0], Ga[1], Ga[2],
                           L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], 2*pi/3)

            dimTag_4 = factory.copy(dimTag_2)
            dimTag_6 = factory.copy(dimTag_2)
            factory.rotate(dimTag_4, Ga[0], Ga[1], Ga[2],
                           L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], -2*pi/3)
            factory.rotate(dimTag_6, Ga[0], Ga[1], Ga[2],
                           L[0]-Ga[0], L[1]-Ga[1], L[2]-Ga[2], 2*pi/3)

            dimTag_24_1 = factory.fuse(
                dimTag, dimTag_6, removeObject=True, removeTool=True)[0]
            dimTag_24_2 = factory.fuse(
                dimTag_4, dimTag_5, removeObject=True, removeTool=True)[0]
            dimTag_24_3 = factory.fuse(
                dimTag_2, dimTag_3, removeObject=True, removeTool=True)[0]

            if mesh_type == "8":
                computeIntersectingVolumes()
                setSize(dimTag_24_1[0][1], h_fine)
                setSize(dimTag_24_2[0][1], h_fine)
                setSize(dimTag_24_3[0][1], h_fine)

            if mesh_type == "2" or mesh_type == "1":
                # dimTag_24_1 is near the x axis.
                # dimTag_24_2 is near the z axis.
                # dimTag_24_2 is near the y axis.
                #
                # The rotation occurs along the z axis.
                #
                # We want to copy and rotate 4 time each.
                # Then we need to fuse the right ones.

                # dimTag_24_1
                dimTag_24_1_1 = factory.copy(dimTag_24_1)
                dimTag_24_1_2 = factory.copy(dimTag_24_1)
                dimTag_24_1_3 = factory.copy(dimTag_24_1)
                factory.rotate(dimTag_24_1_1, 0, 0, 0, 0, 0, 1, pi/2)
                factory.rotate(dimTag_24_1_2, 0, 0, 0, 0, 0, 1, 2*pi/2)
                factory.rotate(dimTag_24_1_3, 0, 0, 0, 0, 0, 1, 3*pi/2)

                # dimTag_24_2
                dimTag_24_2_1 = factory.copy(dimTag_24_2)
                dimTag_24_2_2 = factory.copy(dimTag_24_2)
                dimTag_24_2_3 = factory.copy(dimTag_24_2)
                factory.rotate(dimTag_24_2_1, 0, 0, 0, 0, 0, 1, pi/2)
                factory.rotate(dimTag_24_2_2, 0, 0, 0, 0, 0, 1, 2*pi/2)
                factory.rotate(dimTag_24_2_3, 0, 0, 0, 0, 0, 1, 3*pi/2)

                # dimTag_24_3
                dimTag_24_3_1 = factory.copy(dimTag_24_3)
                dimTag_24_3_2 = factory.copy(dimTag_24_3)
                dimTag_24_3_3 = factory.copy(dimTag_24_3)
                factory.rotate(dimTag_24_3_1, 0, 0, 0, 0, 0, 1, pi/2)
                factory.rotate(dimTag_24_3_2, 0, 0, 0, 0, 0, 1, 2*pi/2)
                factory.rotate(dimTag_24_3_3, 0, 0, 0, 0, 0, 1, 3*pi/2)

                dimTag_4_xmax = factory.fuse(
                    dimTag_24_1,   dimTag_24_3_1, removeObject=True, removeTool=True)[0]
                dimTag_4_xmin = factory.fuse(
                    dimTag_24_1_2, dimTag_24_3_3, removeObject=True, removeTool=True)[0]
                dimTag_4_ymax = factory.fuse(
                    dimTag_24_1_1, dimTag_24_3_2, removeObject=True, removeTool=True)[0]
                dimTag_4_ymin = factory.fuse(
                    dimTag_24_1_3, dimTag_24_3, removeObject=True, removeTool=True)[0]
                dimTag_4_zmax = factory.fuse(
                    dimTag_24_2_1 + dimTag_24_2_2 + dimTag_24_2_3, dimTag_24_2, removeObject=True, removeTool=True)[0]

                if mesh_type == "2":
                    computeIntersectingVolumes()
                    setSize(dimTag_4_xmin[0][1], h_fine)
                    setSize(dimTag_4_xmax[0][1], h_fine)
                    setSize(dimTag_4_ymax[0][1], h_fine)
                    setSize(dimTag_4_ymin[0][1], h_fine)
                    setSize(dimTag_4_zmax[0][1], h_fine)

                if mesh_type == "1":
                    dimTag_4_zmin = factory.copy(dimTag_4_zmax)
                    dimTag_4_xmax_2 = factory.copy(dimTag_4_xmax)
                    dimTag_4_xmin_2 = factory.copy(dimTag_4_xmin)
                    dimTag_4_ymax_2 = factory.copy(dimTag_4_ymax)
                    dimTag_4_ymin_2 = factory.copy(dimTag_4_ymin)

                    # Mirror by the plane of equation "z=0"
                    factory.mirror(dimTag_4_zmin, 0, 0, 1, 0)
                    factory.mirror(dimTag_4_xmax_2, 0, 0, 1, 0)
                    factory.mirror(dimTag_4_xmin_2, 0, 0, 1, 0)
                    factory.mirror(dimTag_4_ymax_2, 0, 0, 1, 0)
                    factory.mirror(dimTag_4_ymin_2, 0, 0, 1, 0)

                    # Fuse the x/y min/max
                    dimTag_1_0 = factory.fuse(
                        dimTag_4_ymax, dimTag_4_ymax_2, removeObject=True, removeTool=True)[0]
                    dimTag_1_1 = factory.fuse(
                        dimTag_4_ymin, dimTag_4_ymin_2, removeObject=True, removeTool=True)[0]
                    dimTag_1_2 = factory.fuse(
                        dimTag_4_xmax, dimTag_4_xmax_2, removeObject=True, removeTool=True)[0]
                    dimTag_1_3 = factory.fuse(
                        dimTag_4_xmin, dimTag_4_xmin_2, removeObject=True, removeTool=True)[0]

                    computeIntersectingVolumes()

                    setSize(dimTag_1_0[0][1], h_fine)
                    setSize(dimTag_1_1[0][1], h_fine)
                    setSize(dimTag_1_2[0][1], h_fine)
                    setSize(dimTag_1_3[0][1], h_fine)

                    setSize(dimTag_4_zmax[0][1], h_fine)
                    setSize(dimTag_4_zmin[0][1], h_fine)


def initialize_gmsh_model(lengthMin, lengthMax, scalingFactor):
    """
    This function initialize the gmsh model. It handles all the technical stuff
    for you: 
    - lengthMin and lengthMax that define how finer the mesh is.
    - algo, ie the type of algorithm used for meshing (algo = 6 by default)
    - ScalingFactor is the coefficient applied to all points corrdinates (4*pi/a by defaults).

    input: 
    lengthMin       (double): Minimum element mesh size.
    lengthMax       (double): Maximum element mesh size.
    scalingFactor   (double):

    output:
    (None) Information is stored in gmsh.model .
    """
    algo = 6  # default
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.ScalingFactor", scalingFactor)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lengthMin)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lengthMax)
    gmsh.option.setNumber("Mesh.SaveElementTagType", 2)
    gmsh.option.setNumber("Mesh.Algorithm", algo)


def build_1_over_48_of_the_BZ(BZ_points):
    """
    This function build 1 over 48 of the brillouin zone. Basically the volume
    delimited  by the usual 1D path inside the BZ.

    input:
    BZ_points       dict: All keys are Ga, L, U, X, W, K. 
                          And each key is np array of 3 coordinates (Kx, Ky, Kz).

    output:         
    None (stored in gmsh.model)
    """

    Ga = BZ_points["Ga"]
    L = BZ_points["L"]
    U = BZ_points["U"]
    X = BZ_points["X"]
    W = BZ_points["W"]
    K = BZ_points["K"]

    # Adding points. ORDER MATTER !
    all_points = np.array([Ga, L, U, X, W, K])
    tag_point = []

    for p in all_points:
        tag_point.append(factory.addPoint(p[0], p[1], p[2]))

    # Adding line. ORDER MATTER !
    tag_line = []

    tag_line.append(factory.addLine(tag_point[0], tag_point[1]))  # Ga -> L
    tag_line.append(factory.addLine(tag_point[1], tag_point[2]))  # L  -> U
    tag_line.append(factory.addLine(tag_point[2], tag_point[3]))  # U  -> X
    tag_line.append(factory.addLine(tag_point[3], tag_point[0]))  # X  -> Ga
    tag_line.append(factory.addLine(tag_point[3], tag_point[4]))  # X  -> W
    tag_line.append(factory.addLine(tag_point[2], tag_point[4]))  # U  -> W
    tag_line.append(factory.addLine(tag_point[4], tag_point[5]))  # W  -> K
    tag_line.append(factory.addLine(tag_point[5], tag_point[1]))  # K  -> L
    tag_line.append(factory.addLine(tag_point[5], tag_point[0]))  # K  -> Ga

    # Adding curve loop.
    tag_curveloop = []

    tag_curveloop.append(factory.addCurveLoop(
        [tag_line[0], tag_line[1], tag_line[2], tag_line[3]]))   # Ga -> L -> U -> X -> Ga
    tag_curveloop.append(factory.addCurveLoop(
        [tag_line[0], -tag_line[7], tag_line[8]]))              # Ga -> L -> K -> Ga
    tag_curveloop.append(factory.addCurveLoop(
        [tag_line[4], -tag_line[5], tag_line[2]]))              # X -> W -> U -> X
    tag_curveloop.append(factory.addCurveLoop(
        [-tag_line[1], -tag_line[7], -tag_line[6], -tag_line[5]]))  # U -> L -> K -> W -> U
    tag_curveloop.append(factory.addCurveLoop(
        [-tag_line[8], -tag_line[6], -tag_line[4], tag_line[3]]))  # Ga -> K -> W -> X -> Ga

    # Adding Surface
    tag_surface = []
    for l in tag_curveloop:
        tag_surface.append(factory.addPlaneSurface([l]))

    # Adding Surface loop
    tag_surfaceloop = factory.addSurfaceLoop(tag_surface)

    # Adding Volume
    Tag = factory.addVolume([tag_surfaceloop])
    factory.synchronize()
    return [(3, Tag)]


def writeMeshFile(filename, view=False, generate=True, dim=3, extension=".msh"):
    """
    Generic function to:
    Either view the geometric model in gmsh gui. (view=True, generate=False)
    Either view the geometric model and the mesh in gmsh gui. (view=True, generate=True)
    Either write the mesh file. (default values)

    input:
    filename    (str):  name of mesh to write or view. Arg extension is always concatenated.
    view        (bool): If you wish to view the model or mesh in the gmsh gui.
    generate    (bool): If you wish to generate the mesh.
    dim         (bool): Dimension of the mesh, either 1, 2, or 3. For meshing volumes keep it to 3.
    extension   (str):  Type of meshfile written. See the docs of gmsh for all available format.
                        extension is always concatenated to filename.
    """
    if(generate):
        model.mesh.generate(dim)
    if(view):
        gmsh.fltk.run()
    else:
        gmsh.write(filename + extension)
    gmsh.finalize()


def computeIntersectingVolumes():
    """
    Compute the intersection between volumes inside the 
    gmsh model (gmsh.model.occ).

    Basically, this function just does a fragment on all volumes.
    """
    factory.synchronize()
    allv = model.getEntities(3)
    ov, ovv = factory.fragment(allv, [])
    factory.synchronize()


def setSize(Tag_refinement, size):
    """
    This function refine the mesh in the object of dimension 3 and of 
    tag equal to "Tag_refinement", with the size "size".
    bounding box defined by bounding_box, with the size "size".

    input:
    Tag_refinement   int: tag of the volume to set finer mesh on.
    size             double: Desired size of the mesh elements.

    output:     
    None (information is stored in gmsh.model)
    """
    xmin, ymin, zmin, xmax, ymax, zmax = factory.getBoundingBox(
        3, Tag_refinement)

    box = np.array(model.getEntitiesInBoundingBox(xmin,
                                                  ymin,
                                                  zmin,
                                                  xmax,
                                                  ymax,
                                                  zmax,
                                                  dim=0))
    model.mesh.setSize(box, size)


if __name__ == "__main__":
    main_BZ_generating_mesh()
