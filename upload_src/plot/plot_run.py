import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Pmw.Pmw_2_0_1.doc.ExampleDemo import title
from matplotlib.pyplot import ylabel

# os.chdir('/home/ding/桌面/work/project_protein/outdata/md__resdata/run_md')





command_str = (f'''
conda activate md
conda activate py312

# 1. 轨迹处理 去周期性 居中处理
# 分析蛋白质-溶剂相互作用; 计算溶剂可及表面积(SASA); 研究水合层结构; 分析离子分布; 保持溶剂环境完整性
echo "0" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -ur compact
echo -e "1\n0" | gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_center.xtc -center -pbc mol -ur compact

# 旋转+平移拟合 (-fit rot+trans) 保持蛋白质内部构象变化不变
# 计算蛋白质内部运动（如RMSD、RMSF）; 分析蛋白质构象变化; 准备结构叠合可视化; 消除整体运动后研究内部动力学
echo -e "1\n1" | gmx trjconv -s md.tpr -f md.xtc -o md_fitted.xtc -fit rot+trans

# 2. 基础分析
echo -e "Backbone\nBackbone" | gmx rms -s md.tpr -f md_center.xtc -o rmsd.xvg
echo -e "C-alpha" | gmx rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg -res
echo -e "Protein" | gmx gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg # 回旋半径(Rg)计算
echo -e "Potential Temperature Pressure Density" | gmx energy -f md.edr -o energy.xvg 
echo -e "Protein" | gmx sasa -s md.tpr -f md_center.xtc -o sasa.xvg -tu ns # 溶剂可及表面积(SASA)计算
'''
               )

#
# import subprocess
# result = subprocess.run(
#     command_str,
#     cwd='./',
#     check=True,
#     shell=True,
#     text=True,
#     encoding='utf-8'
# )

os.chdir('/home/ding/new_disk/md__resdata/run_md/run')
# os.chdir('/home/ding/new_disk/md__resdata/run_rand_md/run')
os.makedirs('plots', exist_ok=True)

def plot_xvg(file_path, output_format='png', dpi=300, line_color='blue',
             line_width=1.5, grid_alpha=0.3, font_size=12, xlabel = "Time (ns)", ylabel = "Value", title = "MD Analysis",
             x_axis_label_reset = 100, show_point = 1000):
    import re
    import math
    """
    绘制 GROMACS 生成的 .xvg 文件
    参数:
        file_path: .xvg 文件路径
        output_format: 输出图片格式 (png, pdf, svg 等)
        dpi: 图片分辨率
        line_color: 线条颜色
        line_width: 线条宽度
        grid_alpha: 网格透明度
        font_size: 字体大小
    """
    # 读取文件内容
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # 提取元数据 (标题, 轴标签)
    # title = "MD Analysis"
    # xlabel = "Time (ns)"
    # ylabel = "Value"

    for line in lines:
        if line.startswith('@ title'):
            title = re.search(r'"([^"]*)"', line).group(1)
        elif line.startswith('@ xaxis  label'):
            xlabel = re.search(r'"([^"]*)"', line).group(1)
            # 自动将 ps 转换为 ns
            if "(ps)" in xlabel:
                xlabel = xlabel.replace("(ps)", "(ns)")
        elif line.startswith('@ yaxis  label'):
            ylabel = re.search(r'"([^"]*)"', line).group(1)

    # 提取数据
    data = []
    for line in lines:
        if not line.startswith(('#', '@')):
            try:
                # 跳过空行
                if line.strip():
                    row = [float(x) for x in line.split()]
                    # 自动将 ps 转换为 ns
                    if "(ps)" in xlabel or "Time (ps)" in xlabel:
                        row[0] /= 1000.0
                    data.append(row)
            except ValueError:
                continue

    if not data:
        print(f"错误: 文件 {file_path} 中没有可用的数据")
        return

    data = np.array(data)

    if data.shape[0] > show_point:
        tmp = np.linspace(0,data.shape[0]-1,show_point)
        tmp = [round(x) for x in tmp]
        data = data[tmp, :]


    time = data[:, 0]
    values = data[:, 1]

    # 缩放X
    max_time = max(time)/x_axis_label_reset
    time = [x / max_time for x in time]

    # 创建图形
    plt.figure(figsize=(10, 6))

    # 绘制数据
    plt.plot(time, values, linewidth=line_width, color=line_color)

    # 设置标题和标签
    plt.title(title, fontsize=font_size + 2, fontweight='bold', pad=20)
    plt.xlabel(xlabel, fontsize=font_size)
    plt.ylabel(ylabel, fontsize=font_size)

    # 设置网格
    plt.grid(alpha=grid_alpha)

    # 设置坐标轴范围
    plt.xlim(min(time), max(time))
    # plt.xlim(xlim[0], xlim[1])
    y_min, y_max = min(values), max(values)
    y_range = y_max - y_min
    plt.ylim(y_min - 0.05 * y_range, y_max + 0.05 * y_range)

    # 设置刻度
    plt.tick_params(axis='both', which='major', labelsize=font_size - 1)

    # 添加紧凑布局
    plt.tight_layout()

    # 生成输出文件名
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_file = f"{base_name}.{output_format}"

    # 保存图片
    if output_format == 'png':
        plt.savefig(f'plots/{output_file}', dpi=dpi, bbox_inches='tight')
    elif output_format == 'pdf':
        plt.savefig(f'plots/{output_file}', bbox_inches='tight')
    print(f"图片已保存为: {output_file}")

    # 显示图片
    plt.show()

    # 关闭图形
    plt.close()

# 使用示例
plot_xvg('rmsd.xvg', 'pdf', ylabel='RMSD(nm)', title = 'RMSD')
plot_xvg('rmsf.xvg', 'pdf', ylabel='RMSF(nm)', title = 'RMSF')
plot_xvg('gyrate.xvg', 'pdf', ylabel='Radius (nm)', title = 'Gyrate') # RG
plot_xvg('energy.xvg', 'pdf', ylabel='(kg/m^3)', title = 'Energy')






import numpy as np
import matplotlib.pyplot as plt

# 读取数据
data = np.loadtxt('sasa.xvg', comments=['#', '@'])
time = data[:, 0]  # 时间 (ns)
sasa = data[:, 1]  # SASA (nm²)

# 创建图形
plt.figure(figsize=(10, 6))
plt.plot(time, sasa, 'g-', linewidth=2)
plt.title('Solvent Accessible Surface Area (SASA)', fontsize=14, pad=20)
plt.xlabel('Time (ns)', fontsize=12)
plt.ylabel('SASA (nm²)', fontsize=12)
plt.grid(alpha=0.3)

# 添加统计信息和关键值
mean_sasa = np.mean(sasa)
min_sasa = np.min(sasa)
max_sasa = np.max(sasa)

plt.axhline(y=mean_sasa, color='r', linestyle='--', alpha=0.7,
            label=f'Mean: {mean_sasa:.1f} nm²')
plt.axhline(y=min_sasa, color='purple', linestyle=':', alpha=0.5,
            label=f'Min: {min_sasa:.1f} nm²')
plt.axhline(y=max_sasa, color='orange', linestyle=':', alpha=0.5,
            label=f'Max: {max_sasa:.1f} nm²')

plt.legend(loc='best')

# 保存并显示
plt.tight_layout()
plt.savefig('plots/sasa_plot.pdf')
plt.show()



# gib 自由能景观
'''
# 常用反应坐标组合
reaction_coord_pairs = [
    ('Rg', 'RMSD'),  # 结构紧凑度 vs 结构偏差
    ('PCA1', 'PCA2'), # 主成分分析的前两个分量
    ('Distance', 'Angle'), # 特定原子间距离和角度
    ('SASA', 'Rg'),   # 溶剂可及表面积 vs 回旋半径
    ('HBonds', 'RMSD') # 氢键数量 vs 结构偏差
]
'''




# 读取文件内容
with open('rmsd.xvg', 'r') as f:
    lines = f.readlines()
# 提取数据
data_rmsd = []
for line in lines:
    if not line.startswith(('#', '@')):
        try:
            # 跳过空行
            if line.strip():
                row = [float(x) for x in line.split()]
                # 自动将 ps 转换为 ns
                data_rmsd.append(row)
        except ValueError:
            continue

with open('gyrate.xvg', 'r') as f:
    lines = f.readlines()
data_gyrate = []
for line in lines:
    if not line.startswith(('#', '@')):
        try:
            # 跳过空行
            if line.strip():
                row = [float(x) for x in line.split()]
                # 自动将 ps 转换为 ns
                data_gyrate.append(row)
        except ValueError:
            continue
data_rmsd = pd.DataFrame(data_rmsd).iloc[:,0:2]
data_gyrate = pd.DataFrame(data_gyrate).iloc[:,0:2]
data_rmsd.columns = ['time', 'rmsd']
data_gyrate.columns = ['time', 'gyrate']
combined = pd.merge(data_rmsd, data_gyrate, on='time', how='inner')
combined.to_csv('output.xvg', sep=' ', index=False, header=False)

'''
gmx sham -f output.xvg -ls gibbs.xpm -nlevels 100 -ngrid 200 -tsham 310 -lsh enthalpy.xpm -lss entropy.xpm # -minfree 0 -maxfree 14
gmx xpm2ps -f gibbs.xpm -o ./plots/gibbs.eps -rainbow red



find . -type f ! -name '#*' ! -size +1G -print0 | \
  tar -czvf backup.tar.gz --null -T -
'''




'''
echo -e "4\n4" | gmx rms -s md.tpr -f md_noPBC.xtc -o protein_rmsd.xvg
echo -e "13\n13" | gmx rms -s md.tpr -f md_noPBC.xtc -o ligand_rmsd.xvg

'''


def plot_multiple_xvg(file_paths, labels, colors, output_format='png', dpi=300,
                      line_width=1.5, grid_alpha=0.3, font_size=12,
                      xlabel="Time (ns)", ylabel="RMSD (nm)", title="RMSD Comparison",
                      show_point=1000, figsize=(10, 6), legend_loc='best', x_axis_label_reset = 100):
    """
    绘制多个 GROMACS 生成的 .xvg 文件在同一张图中

    参数:
        file_paths: .xvg 文件路径列表
        labels: 每条曲线的标签列表
        colors: 每条曲线的颜色列表
        output_format: 输出图片格式 (png, pdf, svg 等)
        dpi: 图片分辨率
        line_width: 线条宽度
        grid_alpha: 网格透明度
        font_size: 字体大小
        xlabel: X轴标签
        ylabel: Y轴标签
        title: 图表标题
        show_point: 最大显示点数 (用于下采样)
        figsize: 图表尺寸
        legend_loc: 图例位置
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import re

    # 确保输入参数长度匹配
    if not (len(file_paths) == len(labels) == len(colors)):
        raise ValueError("file_paths, labels, and colors must have the same length")

    plt.figure(figsize=figsize)
    all_times = []
    all_values = []
    print(1)
    for i, file_path in enumerate(file_paths):
        # 读取文件内容
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # 提取数据
        data = []
        for line in lines:
            if not line.startswith(('#', '@')):
                try:
                    if line.strip():
                        row = [float(x) for x in line.split()]
                        data.append(row)
                except ValueError:
                    continue

        if not data:
            print(f"警告: 文件 {file_path} 中没有可用的数据，跳过")
            continue

        data = np.array(data)
        # 下采样处理
        if data.shape[0] > show_point:
            step = int(data.shape[0] / show_point)
            data = data[::step, :]

        time = data[:, 0]
        values = data[:, 1]
        # 缩放X
        max_time = max(time) / x_axis_label_reset
        time = [x / max_time for x in time]
        # 存储所有数据用于设置坐标轴范围
        all_times.append(time)
        all_values.append(values)

        # 绘制曲线
        plt.plot(time, values, label=labels[i], color=colors[i],
                 linewidth=line_width, alpha=0.8)

    # 设置标题和标签
    plt.title(title, fontsize=font_size + 2, fontweight='bold', pad=20)
    plt.xlabel(xlabel, fontsize=font_size)
    plt.ylabel(ylabel, fontsize=font_size)

    # 设置网格
    plt.grid(alpha=grid_alpha, linestyle='--')

    # 设置坐标轴范围
    if all_times:
        min_time = min([min(t) for t in all_times])
        max_time = max([max(t) for t in all_times])
        plt.xlim(min_time, max_time)

        min_val = min([min(v) for v in all_values])
        max_val = max([max(v) for v in all_values])
        val_range = max_val - min_val
        plt.ylim(max(0, min_val - 0.1 * val_range), max_val + 0.1 * val_range)

    # 设置刻度
    plt.tick_params(axis='both', which='major', labelsize=font_size - 1)

    # 添加图例
    plt.legend(fontsize=font_size, loc=legend_loc, framealpha=0.8)

    # 添加紧凑布局
    plt.tight_layout()

    # 生成输出文件名
    output_file = f"combined_rmsd.{output_format}"

    # 确保plots目录存在
    os.makedirs('plots', exist_ok=True)

    # 保存图片
    save_path = os.path.join('plots', output_file)
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    print(f"图片已保存为: {save_path}")

    # 显示图片
    plt.show()

    # 关闭图形
    plt.close()
file_paths = ["protein_rmsd.xvg", "ligand_rmsd.xvg"]
labels = ["Protein", "Ligand"]
colors = ["#1f77b4", "#ff7f0e"]  # 蓝色和橙色
plot_multiple_xvg(
    file_paths=file_paths,
    labels=labels,
    colors=colors,
    title="Protein and Ligand RMSD Comparison",
    ylabel="RMSD (nm)",
    xlabel="Time (ns)",
    line_width=2.0,
    font_size=12,
    dpi=300,
    show_point=1000,
    legend_loc='upper right'
)
plot_xvg('protein_rmsd.xvg', 'pdf', ylabel='RMSD(nm)', title = 'RMSD')
plot_xvg('ligand_rmsd.xvg', 'pdf', ylabel='RMSD(nm)', title = 'RMSD')