

# 选择并显示结合位点
select binding_site, byres(ligand around 5)

show sticks, binding_site
color palegreen, binding_site



# 标注残基
# label binding_site and name ca, "%s%s" % (resn, resi)
set label_size, 20

# 渲染图像
zoom ligand
ray 2400,2400
png binding_site_view.png, dpi=300

# 显示相互作用
distance hbonds, ligand, binding_site, 3.2, mode=2
set dash_color, purple, hbonds




# 展示配体周围氨基酸 1

# 首先选择配体（假设它是唯一的有机小分子）
select ligand, organic

# 然后选择配体周围5Å内的残基
select binding_site, byres (ligand around 5)

# 显示结果
show sticks, binding_site
show surface, minimized
label binding_site and name ca, resn
color green, binding_site






# 展示配体周围氨基酸 2

# 确保配体可见
show sticks, organic
color yellow, organic

# 选择结合位点（使用方法1）
select ligand, organic
select binding_site, byres (ligand around 5)

# 显示和标注
show sticks, binding_site
show surface, minimized
set surface_color, lightblue, minimized
label binding_site and name ca, resn
color green, binding_site

# 可选：测量距离
distance hbonds, ligand, binding_site, 3.2, mode=2





select binding_residues, byres (ligand around 5)
create binding_site_only, binding_residues and not solvent and not inorganic
save binding_site_residues.pdb, binding_site_only








#################
load protein.pdb, minimized
load ligand.pdbqt, ligand

create binding_site_obj, byres(ligand around 5)

# 显示
as sticks, binding_site_obj
as sticks, ligand
as surface, minimized

# 颜色
color palegreen, binding_site_obj
#color yellow, ligand
set surface_color, lightblue, minimized

# 标签
set label_size, 20
label binding_site_obj and name ca, "%s %s" % (resn, resi)



# 设置距离线样式
set dash_radius, 0.3
set dash_gap, 0.3
set dash_color, red
# # 盐桥
# distance electrostatic, ligand, binding_site_obj, 5.0, 2
#
# # 氢键
# distance h_bonds, ligand, binding_site_obj, 3.5, 1
# color green, h_bonds
#
# # 疏水接触
# distance hydrophobic, ligand, binding_site_obj, 4.5, 1
# color gray, hydrophobic


# set label_screen_point, 1
# set label_position, [0.02, 0.95]  # 左上角
# label binding_site_obj and name ca, "%s %s" % (resn, resi)

# 视角和渲染
zoom ligand, 2
# orient



# 来一个2埃内的所有aa. 并将文本展示于最前
# 创建2Å内的结合口袋残基
create binding_pocket_2A, byres(ligand around 2)

# 显示样式
show sticks, binding_pocket_2A
color cyan, binding_pocket_2A
set stick_radius, 0.15, binding_pocket_2A

select ca_atoms, binding_pocket_2A and name CA
label ca_atoms, "%s %s %s" % (resn, resi, chain)
set label_color, white
set label_size, 18
set label_font_id, 7


save session.pse

