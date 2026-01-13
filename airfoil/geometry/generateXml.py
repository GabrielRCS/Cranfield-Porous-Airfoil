N = 150
ly = 8
lx = 8

maxT = 20.1
vtkSave = 20
imSave = 0.5

for i in range(1, 29):  # 1 through 28
    content = f"""<geometry>
    <location>./geometry/geometry_{i}.dat</location>
    <ResultFile>{i}</ResultFile>
    <N>{N}</N>
    <lx>{lx}</lx>
    <ly>{ly}</ly>
</geometry>
<simulation>
    <maxT>{maxT}</maxT>
    <vtkSave>{vtkSave}</vtkSave>
    <imSave>{imSave}</imSave>
</simulation>
"""

    filename = f"parameters_{i}.xml"
    with open(filename, "w") as f:
        f.write(content)

    print(f"Generated {filename}")
