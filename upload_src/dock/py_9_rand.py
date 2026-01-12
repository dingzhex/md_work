

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from keras.src.utils.file_utils import makedirs
from mpl_toolkits.mplot3d import Axes3D
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# ======================
# 1. 结果提取与整合
# ======================
raw_path = os.getcwd()
os.makedirs('summary_rand', exist_ok=True)
os.chdir('summary_rand')

def parse_vina_results(result_file):
    """解析Vina输出文件，提取结合能和构象信息"""
    # result_file = '/home/ding/桌面/work/project_protein/outdata/aa_site_position/existed/protein_site_GLU17/vina_res/ligand_3d_initial_137_output_file'
    results = []
    current_energy = None
    current_pose = []

    with open(result_file, 'r') as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                if current_energy is not None:
                    results.append({
                        "energy": current_energy,
                        "pose": "".join(current_pose)
                    })
                current_energy = float(line.split()[3])
                current_pose = []
            elif line.startswith("MODEL") or line.startswith("ENDMDL"):
                continue
            else:
                current_pose.append(line)

    # 添加最后一个构象
    if current_energy is not None:
        results.append({
            "energy": current_energy,
            "pose": "".join(current_pose)
        })

    return results


def collect_all_results(results_dir):
    # results_dir = '/home/ding/桌面/work/project_protein/outdata/aa_site_position/existed/protein_site_GLU17'
    """收集所有对接结果"""
    all_results = []
    result_files = [results_dir]

    for file_path in result_files:

        # 解析结果
        docking_results = parse_vina_results(file_path)

        # 记录最佳构象
        best_result = min(docking_results, key=lambda x: x["energy"])


        filename = os.path.basename(file_path)
        parts = filename.split("_")
        ligand = '_'.join(parts[0:4])
        site = os.path.basename(results_dir)

        # 解析结果
        docking_results = parse_vina_results(file_path)
        # 添加到总结果
        all_results.append({
            "ligand": ligand,
            "site": site,
            "file": file_path,
            "best_energy": best_result["energy"],
            "mean_energy": np.mean([r["energy"] for r in docking_results]),
            "std_energy": np.std([r["energy"] for r in docking_results]),
            "num_poses": len(docking_results),
            "best_pose": best_result["pose"]
        })

    return pd.DataFrame(all_results)


# ======================
# 2. 结合能可视化分析
# ======================

def plot_energy_comparison(df):
    """绘制结合能比较图"""
    plt.figure(figsize=(12, 8))

    # 创建分组条形图
    ax = sns.barplot(
        x="ligand",
        y="best_energy",
        hue="site",
        data=df,
        palette="viridis",
        color='.2',
        capsize=.1
    )

    plt.title("Binding Energy Comparison by Site and Ligand", fontsize=16)
    plt.ylabel("Binding Energy (kcal/mol)", fontsize=14)
    plt.xlabel("Ligand", fontsize=14)
    plt.xticks(rotation=45)
    plt.legend(title="Active Site", loc="upper right")

    # 添加数据标签
    for p in ax.patches:
        ax.annotate(f"{p.get_height():.2f}",
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center',
                    xytext=(0, 10),
                    textcoords='offset points',
                    fontsize=9)

    plt.tight_layout()
    plt.savefig("binding_energy_comparison.png")
    plt.close()


# ======================
# 3. 构象稳定性分析
# ======================

def plot_energy_distribution(df):
    """绘制能量分布雷达图"""
    sites = df["site"].unique()
    ligands = df["ligand"].unique()

    # 设置雷达图角度
    angles = np.linspace(0, 2 * np.pi, len(sites), endpoint=False).tolist()
    angles += angles[:1]  # 闭合图形

    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, polar=True)

    for ligand in ligands:
        ligand_data = df[df["ligand"] == ligand]
        values = ligand_data["best_energy"].tolist()
        values += values[:1]  # 闭合图形

        ax.plot(angles, values, linewidth=2, label=ligand)
        ax.fill(angles, values, alpha=0.1)

    # 添加标签
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.degrees(angles[:-1]), sites)

    # 设置径向轴
    min_energy = df["best_energy"].min() - 1
    max_energy = df["best_energy"].max() + 1
    ax.set_ylim(min_energy, max_energy)
    ax.set_rlabel_position(0)

    plt.title("Binding Energy Distribution Across Sites", fontsize=16)
    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    plt.tight_layout()
    plt.savefig("energy_distribution_radar.png")
    plt.close()


# ======================
# 4. 聚类分析
# ======================

def cluster_analysis(df):
    """对结果进行聚类分析"""
    # 准备数据矩阵：配体×位点
    # 创建数据矩阵
    matrix = df.pivot(index="ligand", columns="site", values="best_energy")

    # 标准化数据
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)

    # 执行层次聚类
    row_linkage = linkage(scaled_data, method='ward', optimal_ordering=True)
    col_linkage = linkage(scaled_data.T, method='ward', optimal_ordering=True)

    # 创建图形和轴 - 使用GridSpec精确控制布局
    fig = plt.figure(figsize=(15, 12))

    # 使用GridSpec定义布局
    gs = plt.GridSpec(4, 4, figure=fig,
                      height_ratios=[0.1, 0.8, 0.1, 0.1],
                      width_ratios=[0.8, 0.1, 0.1, 0.1],
                      hspace=0.05, wspace=0.05)

    # 行树状图 (顶部)
    ax_row = fig.add_subplot(gs[0, 0])
    dendrogram(
        row_linkage,
        orientation='top',
        labels=matrix.index.tolist(),
        color_threshold=0,
        above_threshold_color='k'
    )
    ax_row.set_axis_off()

    # 列树状图 (右侧)
    ax_col = fig.add_subplot(gs[1, 1])
    dendrogram(
        col_linkage,
        orientation='right',
        labels=matrix.columns.tolist(),
        color_threshold=0,
        above_threshold_color='k'
    )
    ax_col.set_axis_off()

    # 主热图
    ax_heatmap = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        matrix,
        cmap="viridis_r",
        annot=True,
        fmt=".2f",
        cbar_kws={'label': 'Binding Energy (kcal/mol)'},
        ax=ax_heatmap
    )

    # 添加颜色条 (单独创建轴)
    cbar_ax = fig.add_subplot(gs[1, 2])
    fig.colorbar(ax_heatmap.collections[0], cax=cbar_ax, label='Binding Energy (kcal/mol)')

    # 添加标题
    plt.suptitle("Hierarchical Clustering of Binding Affinities", fontsize=16, y=0.95)

    # 调整布局并保存
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # 为标题留出空间
    plt.savefig("binding_affinity_clustering.png", bbox_inches='tight', dpi=300)
    # plt.show()
    plt.close()


# ======================
# 5. 多维尺度分析 (MDS)
# ======================

def multidimensional_scaling(df):
    """执行多维尺度分析"""
    # 创建配体-位点特征矩阵
    features = df.pivot(index="ligand", columns="site", values="best_energy").fillna(0)

    # 标准化数据
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    # 执行PCA降维
    pca = PCA()
    components = pca.fit_transform(scaled_features)

    # 创建3D散点图
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # 绘制配体点
    for i, ligand in enumerate(features.index):
        ax.scatter(
            components[i, 0],
            components[i, 1],
            components[i, 2],
            s=200,
            label=ligand,
            alpha=0.8
        )
        ax.text(
            components[i, 0],
            components[i, 1],
            components[i, 2],
            ligand,
            fontsize=12
        )

    # 设置轴标签
    ax.set_xlabel('PC1 ', fontsize=12)
    ax.set_ylabel('PC2 ', fontsize=12)
    ax.set_zlabel('PC3 ', fontsize=12)

    plt.title("3D Projection of Ligand-Site Binding Patterns", fontsize=16)
    plt.tight_layout()
    plt.savefig("mds_3d_projection.png")
    plt.close()

    return components


def multidimensional_scaling_2d(df):
    """执行多维尺度分析"""
    # 创建配体-位点特征矩阵
    features = df.pivot(index="ligand", columns="site", values="best_energy").fillna(0)

    # 标准化数据
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    # 执行PCA降维
    pca = PCA()
    components = pca.fit_transform(scaled_features)

    # 创建3D散点图
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)

    # 绘制配体点
    for i, ligand in enumerate(features.index):
        ax.scatter(
            components[i, 0],
            components[i, 1],

            s=200,
            label=ligand,
            alpha=0.8
        )
        ax.text(
            components[i, 0],
            components[i, 1],

            ligand,
            fontsize=12
        )

    # 设置轴标签
    ax.set_xlabel('PC1 ', fontsize=12)
    ax.set_ylabel('PC2 ', fontsize=12)

    plt.title("2D Projection of Ligand-Site Binding Patterns", fontsize=16)
    plt.tight_layout()
    plt.savefig("mds_2d_projection.png")
    plt.close()

    return components



# ======================
# 6. 关键相互作用分析
# ======================

def analyze_key_interactions(df, key_residues):
    """分析关键残基的相互作用"""
    interaction_results = []

    for _, row in df.iterrows():
        pose = row["best_pose"]
        interactions = {res: 0 for res in key_residues}

        # 简化版相互作用分析 - 实际中应使用专业工具
        for line in pose.split("\n"):
            if "ATOM" in line or "HETATM" in line:
                for residue in key_residues:
                    if residue in line:
                        # 检查距离 - 实际中应计算原子距离
                        interactions[residue] += 1

        interaction_results.append({
            "ligand": row["ligand"],
            "site": row["site"],
            **interactions
        })

    return pd.DataFrame(interaction_results)


# ======================
# 7. 综合评分与排序
# ======================

def calculate_composite_score(df, key_residues):
    """计算综合评分"""
    # 1. 归一化结合能
    df["energy_score"] = 1 - (df["best_energy"] - df["best_energy"].min()) / (
            df["best_energy"].max() - df["best_energy"].min())

    # 2. 归一化能量稳定性
    df["stability_score"] = 1 - (df["std_energy"] / df["std_energy"].max())

    # 3. 关键相互作用评分
    interaction_df = analyze_key_interactions(df, key_residues)
    interaction_scores = []

    for _, row in interaction_df.iterrows():
        score = 0
        for res in key_residues:
            # 每个关键残基最高贡献0.2分
            score += min(0.2, row[res] * 0.05)
        interaction_scores.append(score)

    df["interaction_score"] = interaction_scores

    # 4. 综合评分 (权重可调整)
    df["composite_score"] = (
            0.5 * df["energy_score"] +
            0.3 * df["stability_score"] +
            0.2 * df["interaction_score"]
    )

    # 按综合评分排序
    df = df.sort_values("composite_score", ascending=False)

    return df


# ======================
# 主工作流程
# ======================





if __name__ == "__main__":


    # main()
    # 配置参数
    RESULTS_DIR = '/home/ding/桌面/work/project_protein/outdata/aa_site_randsearch/res'
    # KEY_RESIDUES = ["Lys181", "Glu230", "Asp293", "Arg25", "Arg86",
    # "Lys14", "Glu17"]  # 替换为实际关键残基

    # 1. 收集所有结果
    print("收集对接结果...")
    files =  result_files = glob.glob(os.path.join(RESULTS_DIR,'*output*'))

    res_df = []

    for file in files:
        results_df = collect_all_results(file)
        res_df.append(results_df)

    combined_df = pd.concat(res_df, ignore_index=True)
    combined_df.to_csv("all_docking_results.csv", index=False)
    combined_df.iloc[:,:7].to_csv("all_docking_results_simple.csv", index=False)



    os.chdir(raw_path)
    os.makedirs('outdata/md_rand', exist_ok=True)
    os.makedirs('outdata/md_rand/protein_ligand', exist_ok=True)
    combined_df = combined_df.sort_values('best_energy', ascending=True)

    sel_file = combined_df.iloc[0,2]

    import subprocess
    command_str = (f'''
        cp {sel_file} ./ 
        vina_split --input *output*
        '''
    )
    # 使用单字符串命令 + shell=True
    result = subprocess.run(
        command_str,
        cwd='outdata/md_rand/protein_ligand',
        check=True,
        shell=True,
        text=True,
        encoding='utf-8'
    )

    import os
    run_path = 'outdata/md_rand/run'
    os.makedirs(run_path, exist_ok = True)
    all_lig_file = os.listdir('outdata/md_rand/protein_ligand')
    all_lig_file.sort()
    try:
        idx = all_lig_file[1].split('file_ligand_')[1].split('.pdbqt')[0]
        if int(idx) == 1:
            print('成功选中最优构象')

        command_str = (f'''
            cp {all_lig_file[1]} ../run/ligand.pdb
            cp ../../mm_minimized/minimized.pdb ../run/protein.pdb
            '''
                       )
        # 使用单字符串命令 + shell=True
        result = subprocess.run(
            command_str,
            cwd='outdata/md_rand/protein_ligand',
            check=True,
            shell=True,
            text=True,
            encoding='utf-8'
        )

    except Exception as e:
        print(e)



