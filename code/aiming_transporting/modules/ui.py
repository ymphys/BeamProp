import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button
from PyQt5.QtWidgets import QApplication, QFileDialog
import sys

# 初始化QApplication
if not QApplication.instance():
    app = QApplication(sys.argv)

# 设置中文字体和负号显示（如有需要）
matplotlib.rcParams['font.sans-serif'] = ['Heiti TC', 'STHeiti', 'SimHei', 'Arial Unicode MS']
matplotlib.rcParams['axes.unicode_minus'] = False

def get_user_inputs():
    """
    弹出matplotlib界面，获取用户输入参数。
    返回: dict, 包含所有输入参数和模式
    """
    fig, ax = plt.subplots(figsize=(10, 7))  # 增大窗口尺寸
    plt.subplots_adjust(bottom=0.4)
    plt.title("电子束瞄准与注量率估算：参数输入")
    plt.axis('off')

    # 两列布局参数
    left1, left2 = 0.18, 0.70  # 左右列整体右移
    width, height, vgap = 0.28, 0.05, 0.07
    # 将bottoms整体下移，避免顶部被遮挡
    bottoms = [0.62, 0.54, 0.46, 0.38]

    # 左列：能量、磁场、距离、角度
    E_box = TextBox(plt.axes([left1, bottoms[0], width, height]), '电子束能量 [MeV]:', initial="35")
    B_box = TextBox(plt.axes([left1, bottoms[1], width, height]), '磁场 B [g]:', initial="0.25")
    R_box = TextBox(plt.axes([left1, bottoms[2], width, height]), '目标距离 R [km]:', initial="15")
    Psi_box = TextBox(plt.axes([left1, bottoms[3], width, height]), '目标角度 Ψ [rad]:', initial="0.5")

    # 右列：密度、电流、发射度类型、初始半径
    ni_box = TextBox(plt.axes([left2, bottoms[0], width, height]), r'背景离子密度 $n_i$ [cm$^{-3}$]:', initial="1e5")
    Ib_box = TextBox(plt.axes([left2, bottoms[1], width, height]), r'电子束电流 $I_b$ [A]:', initial="1")
    em_type_box = TextBox(plt.axes([left2, bottoms[2], width, height]), '发射度类型("Ti", "C", or "LP"):', initial="Ti")
    R0_box = TextBox(plt.axes([left2, bottoms[3], width, height]), r'初始束流包络半径 $R_0$ [cm]:', initial="5")

    # 批量仿真按钮放在下方中间
    batch_button_ax = plt.axes([0.35, 0.22, 0.3, 0.06])
    batch_button = Button(batch_button_ax, '批量仿真(选择CSV文件)')

    # 确认按钮放在上方
    button_ax = plt.axes([0.4, 0.30, 0.2, 0.06])
    button = Button(button_ax, '确认')

    user_inputs = {}
    done = {'flag': False}

    def submit(event=None):
        try:
            user_inputs['E'] = float(E_box.text)
            if user_inputs['E'] <= 0:
                raise ValueError('电子束能量必须为正')
            user_inputs['B'] = float(B_box.text)
            if user_inputs['B'] < 0:
                raise ValueError('磁场不能为负')
            user_inputs['R'] = float(R_box.text)
            if user_inputs['R'] <= 0:
                raise ValueError('目标距离必须为正')
            user_inputs['Psi'] = float(Psi_box.text)
            if not (0 < user_inputs['Psi'] < 3.2):
                raise ValueError('目标角度应在0到π之间')
            user_inputs['n_i'] = float(ni_box.text)
            if user_inputs['n_i'] < 0:
                raise ValueError('背景离子密度不能为负')
            user_inputs['Ib'] = float(Ib_box.text)
            if user_inputs['Ib'] <= 0:
                raise ValueError('电子束电流必须为正')
            user_inputs['R0'] = float(R0_box.text)
            if user_inputs['R0'] <= 0:
                raise ValueError('电子束初始半径必须为正')
            user_inputs['em_type'] = em_type_box.text.strip()
            if user_inputs['em_type'] not in ('Ti', 'C', 'LP'):
                raise ValueError('发射度类型必须为 "Ti", "C", 或 "LP"')
            user_inputs['mode'] = 'single'
            done['flag'] = True
            plt.close()
        except ValueError as ve:
            print(f"参数错误: {ve}")
        except Exception:
            print("请输入有效的数值参数")

    def select_csv(event=None):
        file_path, _ = QFileDialog.getOpenFileName(
            None,
            "选择CSV文件",
            "",
            "CSV Files (*.csv)"
        )
        if file_path:
            user_inputs['mode'] = 'batch'
            user_inputs['csv_path'] = file_path
            done['flag'] = True
            plt.close()

    button.on_clicked(submit)
    batch_button.on_clicked(select_csv)

    def on_key(event):
        if event.key == 'enter':
            submit()
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.show()
    if done['flag']:
        return user_inputs
    else:
        return None
