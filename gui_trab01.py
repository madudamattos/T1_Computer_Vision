import sys
import matplotlib.pyplot as plt
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QLabel, QWidget, QLineEdit, QHBoxLayout, QVBoxLayout, QPushButton,QGroupBox
from PyQt5.QtGui import QDoubleValidator
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D
from numpy import array
from stl import mesh

###### Crie suas funções de translação, rotação, criação de referenciais, plotagem de setas e qualquer outra função que precisar

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        #definindo as variaveis
        self.set_variables()
        #Ajustando a tela    
        self.setWindowTitle("Grid Layout")
        self.setGeometry(100, 100,1280 , 720)
        self.setup_ui()

    def set_variables(self):
        self.objeto_original = defineObj() 
        self.objeto = self.objeto_original
        self.cam_original = np.array([[1,0,0,0],
                                      [0,1,0,0],
                                      [0,0,1,0], 
                                      [0,0,0,1]])
        T1 = translate(30,-10,2)
        R1 = x_rotate(-90)
        R2 = z_rotate(90)
        self.cam = T1 @ R2 @ R1 @ self.cam_original
        self.cam_quivers = [] # ajustar
        self.px_base = 1280  
        self.px_altura = 720 
        self.dist_foc = 35 
        self.stheta = 0 
        self.ccd = [36,24] 
        self.projection_matrix = np.array([[1,0,0,0], 
                                           [0,1,0,0], 
                                           [0,0,1,0]])
        self.ox = self.px_base/2 
        self.oy = self.px_altura/2 


    def setup_ui(self):
        # Criar o layout de grade
        grid_layout = QGridLayout()

        # Criar os widgets
        line_edit_widget1 = self.create_world_widget("Ref mundo")
        line_edit_widget2  = self.create_cam_widget("Ref camera")
        line_edit_widget3  = self.create_intrinsic_widget("Params instr")

        self.canvas = self.create_matplotlib_canvas()

        # Adicionar os widgets ao layout de grade
        grid_layout.addWidget(line_edit_widget1, 0, 0)
        grid_layout.addWidget(line_edit_widget2, 0, 1)
        grid_layout.addWidget(line_edit_widget3, 0, 2)
        grid_layout.addWidget(self.canvas, 1, 0, 1, 3)

        # Criar um widget para agrupar o botão de reset
        reset_widget = QWidget()
        reset_layout = QHBoxLayout()
        reset_widget.setLayout(reset_layout)

        # Criar o botão de reset vermelho
        reset_button = QPushButton("Reset")
        reset_button.setFixedSize(50, 30)  # Define um tamanho fixo para o botão (largura: 50 pixels, altura: 30 pixels)
        style_sheet = """
            QPushButton {
                color : white ;
                background: rgba(255, 127, 130,128);
                font: inherit;
                border-radius: 5px;
                line-height: 1;
            }
        """
        reset_button.setStyleSheet(style_sheet)
        reset_button.clicked.connect(self.reset_canvas)
        #reset_button.clicked.connect(lambda: self.reset_canvas)


        # Adicionar o botão de reset ao layout
        reset_layout.addWidget(reset_button)

        # Adicionar o widget de reset ao layout de grade
        grid_layout.addWidget(reset_widget, 2, 0, 1, 3)

        # Criar um widget central e definir o layout de grade como seu layout
        central_widget = QWidget()
        central_widget.setLayout(grid_layout)
        
        # Definir o widget central na janela principal
        self.setCentralWidget(central_widget)

    def create_intrinsic_widget(self, title):
        # Criar um widget para agrupar os QLineEdit
        line_edit_widget = QGroupBox(title)
        line_edit_layout = QVBoxLayout()
        line_edit_widget.setLayout(line_edit_layout)

        # Criar um layout de grade para dividir os QLineEdit em 3 colunas
        grid_layout = QGridLayout()

        line_edits = []
        labels = ['n_pixels_base:', 'n_pixels_altura:', 'ccd_x:', 'ccd_y:', 'dist_focal:', 'sθ:']  # Texto a ser exibido antes de cada QLineEdit

        # Adicionar widgets QLineEdit com caixa de texto ao layout de grade
        for i in range(1, 7):
            line_edit = QLineEdit()
            label = QLabel(labels[i-1])
            validator = QDoubleValidator()  # Validador numérico
            line_edit.setValidator(validator)  # Aplicar o validador ao QLineEdit
            grid_layout.addWidget(label, (i-1)//2, 2*((i-1)%2))
            grid_layout.addWidget(line_edit, (i-1)//2, 2*((i-1)%2) + 1)
            line_edits.append(line_edit)

        # Criar o botão de atualização
        update_button = QPushButton("Atualizar")

        ##### Você deverá criar, no espaço reservado ao final, a função self.update_params_intrinsc ou outra que você queira 
        # Conectar a função de atualização aos sinais de clique do botão
        update_button.clicked.connect(lambda: self.update_params_intrinsc(line_edits))

        # Adicionar os widgets ao layout do widget line_edit_widget
        line_edit_layout.addLayout(grid_layout)
        line_edit_layout.addWidget(update_button)

        # Retornar o widget e a lista de caixas de texto
        return line_edit_widget
    
    def create_world_widget(self, title):
        # Criar um widget para agrupar os QLineEdit
        line_edit_widget = QGroupBox(title)
        line_edit_layout = QVBoxLayout()
        line_edit_widget.setLayout(line_edit_layout)

        # Criar um layout de grade para dividir os QLineEdit em 3 colunas
        grid_layout = QGridLayout()
 
        line_edits = []
        labels = ['X(move):', 'X(angle):', 'Y(move):', 'Y(angle):', 'Z(move):', 'Z(angle):']  # Texto a ser exibido antes de cada QLineEdit

        # Adicionar widgets QLineEdit com caixa de texto ao layout de grade
        for i in range(1, 7):
            line_edit = QLineEdit()
            label = QLabel(labels[i-1])
            validator = QDoubleValidator()  # Validador numérico
            line_edit.setValidator(validator)  # Aplicar o validador ao QLineEdit
            grid_layout.addWidget(label, (i-1)//2, 2*((i-1)%2))
            grid_layout.addWidget(line_edit, (i-1)//2, 2*((i-1)%2) + 1)
            line_edits.append(line_edit)

        # Criar o botão de atualização
        update_button = QPushButton("Atualizar")

        ##### Você deverá criar, no espaço reservado ao final, a função self.update_world ou outra que você queira 
        # Conectar a função de atualização aos sinais de clique do botão
        update_button.clicked.connect(lambda: self.update_world(line_edits))
        # Adicionar os widgets ao layout do widget line_edit_widget
        line_edit_layout.addLayout(grid_layout)
        line_edit_layout.addWidget(update_button)

        # Retornar o widget e a lista de caixas de texto
        return line_edit_widget

    def create_cam_widget(self, title):
        # Criar um widget para agrupar os QLineEdit
        line_edit_widget = QGroupBox(title)
        line_edit_layout = QVBoxLayout()
        line_edit_widget.setLayout(line_edit_layout)

        # Criar um layout de grade para dividir os QLineEdit em 3 colunas
        grid_layout = QGridLayout()

        line_edits = []
        labels = ['X(move):', 'X(angle):', 'Y(move):', 'Y(angle):', 'Z(move):', 'Z(angle):']  # Texto a ser exibido antes de cada QLineEdit

        # Adicionar widgets QLineEdit com caixa de texto ao layout de grade
        for i in range(1, 7):
            line_edit = QLineEdit()
            label = QLabel(labels[i-1])
            validator = QDoubleValidator()  # Validador numérico
            line_edit.setValidator(validator)  # Aplicar o validador ao QLineEdit
            grid_layout.addWidget(label, (i-1)//2, 2*((i-1)%2))
            grid_layout.addWidget(line_edit, (i-1)//2, 2*((i-1)%2) + 1)
            line_edits.append(line_edit)

        # Criar o botão de atualização
        update_button = QPushButton("Atualizar")

        ##### Você deverá criar, no espaço reservado ao final, a função self.update_cam ou outra que você queira 
        # Conectar a função de atualização aos sinais de clique do botão
        update_button.clicked.connect(lambda: self.update_cam(line_edits))

        # Adicionar os widgets ao layout do widget line_edit_widget
        line_edit_layout.addLayout(grid_layout)
        line_edit_layout.addWidget(update_button)

        # Retornar o widget e a lista de caixas de texto
        return line_edit_widget

    def create_matplotlib_canvas(self):
        # Criar um widget para exibir os gráficos do Matplotlib
        canvas_widget = QWidget()
        canvas_layout = QHBoxLayout()
        canvas_widget.setLayout(canvas_layout)

        # Criar um objeto FigureCanvas para exibir o gráfico 2D
        self.fig1, self.ax1 = plt.subplots()
        self.ax1.set_title("Imagem")
        self.canvas1 = FigureCanvas(self.fig1)

        ##### Falta acertar os limites do eixo X
        self.ax1.set_xlim([0,self.px_base])
        
        ##### Falta acertar os limites do eixo Y
        self.ax1.set_ylim([self.px_altura,0])

        ##### Você deverá criar a função de projeção 
        object_2d = self.projection_2d()

        ##### Falta plotar o object_2d que retornou da projeção

        self.ax1.grid('True')
        self.ax1.set_aspect('equal')  
        self.ax1.plot(object_2d[0,:],object_2d[1,:])
        
        canvas_layout.addWidget(self.canvas1)

        # Criar um objeto FigureCanvas para exibir o gráfico 3D
        self.fig2 = plt.figure()
        self.ax2 = self.fig2.add_subplot(111, projection='3d')
        
        ##### Falta plotar o seu objeto 3D e os referenciais da câmera e do mundo

        # Plot do referencial do mundo
        self.plotWorld()

        # Plot da câmera 
        self.plotCam()

        # Plot do objeto 3D
        self.plotObj()

        # Configurando o eixo 3D
        #set_view_object_cam(self.ax2)

        # Configurando o eixo 3D
        lower_lim = -30
        upper_lim = 30

        self.ax2.set_xlim([lower_lim,upper_lim])
        self.ax2.set_ylim([lower_lim,upper_lim])
        self.ax2.set_zlim([-10,upper_lim])

        self.canvas2 = FigureCanvas(self.fig2)
        canvas_layout.addWidget(self.canvas2)


        # Retornar o widget de canvas
        return canvas_widget

    ##### Você deverá criar as suas funções aqui
    
    def plotObj(self):
        # Plotando os pontos e desenhando as linhas
        self.ax2.plot(self.objeto[0,:],self.objeto[1,:],self.objeto[2,:],'r')
        return 
    
    def plotWorld(self):
        # base vector values
        base = np.array([
                         [1, 0, 0],  
                         [0, 1, 0],  
                         [0, 0, 1],  
                         [0, 0, 0]   
                    ]) 

        # origin point
        origin = np.array([[0],[0],[0],[1]])

        draw_arrows(origin,base,self.ax2,3)
        return 
    
    def plotCam(self):
        draw_arrows(self.cam[:, -1], self.cam[:, 0:3], self.ax2,3)
        return

    def update_params_intrinsc(self, line_edits):
        try:

            if line_edits[0].text() != '': 
                self.px_base = float(line_edits[0].text())

            if line_edits[1].text() != '': 
                self.px_altura = float(line_edits[1].text())
            
            if line_edits[2].text() != '': 
                self.ccd[0] = float(line_edits[2].text())
            
            if line_edits[3].text() != '':
                self.ccd[1] = float(line_edits[3].text())
            
            if line_edits[4].text() != '':
                self.dist_foc = float(line_edits[4].text())
            
            if line_edits[5].text() != '':
                self.stheta = float(line_edits[5].text())
            
            
            self.update_canvas()
            #atualiza a projeção 

            # Limpa as caixas de texto depois de atualizar o valor
            for edit in line_edits:
                edit.clear()

        except Exception as e:
            print("Erro ao atualizar parametros intrinsecos:", e)

        return 

    def update_world(self,line_edits):
        try:
            # atualiza as translacoes
            tx = convert_value(line_edits[0].text())
            ty = convert_value(line_edits[2].text())
            tz = convert_value(line_edits[4].text())
            T = translate(tx, ty, tz)
            
            # atualiza as rotações
            rx = convert_value(line_edits[1].text())
            ry = convert_value(line_edits[3].text())
            rz = convert_value(line_edits[5].text())

            Rx = x_rotate(rx)
            Ry = y_rotate(ry)
            Rz = z_rotate(rz)
            R = Rz@Ry@Rx
            
            # multiplica as transformações pela camera
            self.cam = T@R@(self.cam)
            # arredondamento para limpar a impressão da matiz, escondendo os erros muito pequenos (e-16)
            #self.cam = np.round(self.cam, decimals=6)
            self.update_canvas() 
        except Exception as e:
            print("Erro ao atualizar matriz mundo->camera:", e)      

        # Limpa as caixas de texto depois de atualizar o valor
        for edit in line_edits:
            edit.clear()     
        return

    def update_cam(self,line_edits):

        tx = convert_value(line_edits[0].text())
        ty = convert_value(line_edits[2].text())
        tz = convert_value(line_edits[4].text())
        self.cam = translate_cam(tx, ty, tz, self.cam)
        
        rx = convert_value(line_edits[1].text())
        self.cam = x_rotate_cam(rx, self.cam)
        ry = convert_value(line_edits[3].text())
        self.cam = y_rotate_cam(ry, self.cam)

        rz = convert_value(line_edits[5].text())
        self.cam = z_rotate_cam(rz, self.cam)

        # Arredondamento para limpar a impressão da matriz
        self.cam = np.round(self.cam, decimals=6)

        self.update_canvas() 

        # Limpa as caixas de texto depois de atualizar o valor
        for edit in line_edits:
            edit.clear()
        
        return 

        
    def projection_2d(self):  
        M_intrinsics = self.generate_intrinsic_params_matrix()
        M_projection = self.projection_matrix
        M_extrinsics = np.linalg.inv(self.cam)
        M = M_intrinsics @ M_projection @ M_extrinsics

        # Coordenadas do objeto no sistema da câmera
        obj_cam = M_extrinsics @ self.objeto

        # Faz uma máscara para verificar se o objeto está visivel
        mask = obj_cam[2, :] > 0
        obj_visible = self.objeto[:, mask]

        if obj_visible.shape[1] == 0:
            return np.zeros((3, 0))

        obj_projection = M @ obj_visible
        obj_projection = obj_projection / obj_projection[2, :]
        return obj_projection
    
    def generate_intrinsic_params_matrix(self):
        # Matriz intrínseca
        Ox = self.px_base/2
        Oy = self.px_base/2

        Sx = self.px_base/self.ccd[0]
        Sy = self.px_altura/self.ccd[1]
        
        f = self.dist_foc
        Oh = self.stheta

        M_intrinsics = np.array([[f*Sx, f*Oh, Ox],
                                 [   0, f*Sy, Oy],
                                 [   0,    0,  1]])
                
        return M_intrinsics
    
    def update_canvas(self):
        # Limpa o que estava anteriormente em ax1
        self.ax1.cla()
        self.ax1.set_title("Imagem")

        object_2d = self.projection_2d()
        self.ax1.set_ylim([self.px_altura,0])
        self.ax1.set_xlim([0,self.px_base])
        self.ax1.plot(object_2d[0,:],object_2d[1,:])
        
        self.fig1.canvas.draw_idle()
        self.fig1.canvas.flush_events()
        
        
        self.ax2.cla()
        self.ax2.set_title("Objeto") 

        self.plotWorld()
        self.plotCam()
        self.plotObj()

        set_view_object_cam(self.ax2)

        self.fig2.canvas.draw_idle()
        self.fig2.canvas.flush_events()
        return 

    def reset_canvas(self):
        # Volta as variáveis para o valor original
        self.set_variables()

        # Atualiza o canva
        self.update_canvas()
        return
    
# funçõe de criação de Objeto 3D
def defineObj():
    # Objeto 3D 
    pokemon = mesh.Mesh.from_file('models/poliwag.stl')

    # Transforma as coordenadas
    x = pokemon.x.flatten()
    y = pokemon.y.flatten()
    z = pokemon.z.flatten()

    # Esse é o objeto 3D em coordenadas homogeneas
    pokemon = np.array([x.T, y.T, z.T, np.ones(x.size)])

    # Ajustar o obj inicial
    M = np.array([[4, 0, 0, 0], [0, 4, 0, 0], [0, 0, 4, 0], [0,0,0,1]])
    
    R = z_rotate(20)

    R = y_rotate(20)@ R
    
    pokemon = R@pokemon

    pokemon = M@pokemon

    return pokemon

# Funçõe pra convertir valor de interface
def convert_value(text):
    if text == "": return float(0)
    return float(text)

# Funções de rotação
def x_rotate(angle_x):
    return np.array([[1,0,0,0],
                        [0, np.cos(np.deg2rad(angle_x)),-np.sin(np.deg2rad(angle_x)),0],
                        [0, np.sin(np.deg2rad(angle_x)), np.cos(np.deg2rad(angle_x)),0],
                        [0,0,0,1]])
def y_rotate(angle_y):
    return np.array([[np.cos(np.deg2rad(angle_y)),0, np.sin(np.deg2rad(angle_y)),0],
                        [0,1,0,0],
                        [-np.sin(np.deg2rad(angle_y)), 0, np.cos(np.deg2rad(angle_y)),0],
                        [0,0,0,1]])
def z_rotate(angle_z):
    return np.array([[np.cos(np.deg2rad(angle_z)),-np.sin(np.deg2rad(angle_z)),0,0],
                        [np.sin(np.deg2rad(angle_z)),np.cos(np.deg2rad(angle_z)),0,0],
                        [0,0,1,0],
                        [0,0,0,1]])

def x_rotate_cam(xangle, M_x):
    M_inv = np.linalg.inv(M_x)
    cam_orig = np.dot(M_inv, M_x)
    Rx = np.array([[1,0,0,0],
                [0, np.cos(np.deg2rad(xangle)),-np.sin(np.deg2rad(xangle)),0],
                [0, np.sin(np.deg2rad(xangle)), np.cos(np.deg2rad(xangle)),0],
                [0,0,0,1]])
    M_RxO = Rx@cam_orig
    M_Rx = M_x@M_RxO
    return M_Rx
    
def y_rotate_cam(yangle, M_y):
    M_inv = np.linalg.inv(M_y)
    cam_orig = np.dot(M_inv, M_y)
    Ry = np.array([[np.cos(np.deg2rad(yangle)),0, np.sin(np.deg2rad(yangle)),0],
                    [0,1,0,0],
                    [-np.sin(np.deg2rad(yangle)), 0, np.cos(np.deg2rad(yangle)),0],
                    [0,0,0,1]])
    M_RyO = cam_orig@Ry
    M_Ry = M_y@M_RyO
    return M_Ry
    
def z_rotate_cam(zangle, M_z):
    M_inv = np.linalg.inv(M_z)
    cam_orig = np.dot(M_inv, M_z)
    Rz = np.array([[np.cos(np.deg2rad(zangle)),-np.sin(np.deg2rad(zangle)),0,0],
                    [np.sin(np.deg2rad(zangle)),np.cos(np.deg2rad(zangle)),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])
    M_RzO = cam_orig@Rz
    M_Rz = M_z@M_RzO
    return M_Rz

# Funções de translação
def translate(x,y,z):
    return np.array([[1,0,0,x],
                        [0,1,0,y],
                        [0,0,1,z],
                        [0,0,0,1]])
    

def translate_cam(x,y,z, M_T):
    M_inv_T = np.linalg.inv(M_T)
    cam_orig = np.dot(M_inv_T, M_T)
    T = np.array([[1,0,0,x],
                    [0,1,0,y],
                    [0,0,1,z],
                    [0,0,0,1]])
    M_TO = cam_orig@T
    M_acc = M_T@M_TO
    return M_acc

# Função emprestada do material da professora
def draw_arrows(point,base,axis,length=0.5):
    # Plot vector of x-axis
    qx = axis.quiver(point[0],point[1],point[2],base[0,0],base[1,0],base[2,0],color='red',pivot='tail',  length=length)
    # Plot vector of y-axis
    qy = axis.quiver(point[0],point[1],point[2],base[0,1],base[1,1],base[2,1],color='green',pivot='tail',  length=length)
    # Plot vector of z-axis
    qz = axis.quiver(point[0],point[1],point[2],base[0,2],base[1,2],base[2,2],color='blue',pivot='tail',  length=length)
    return [qx, qy, qz]

# Função para ajustar o eixo de coordenadas do objeto 3D em relação à mudança de posição ou orientação da câmera
def set_view_object_cam(ax):

    x_lim = ax.get_xlim3d()
    y_lim = ax.get_ylim3d()
    z_lim = ax.get_zlim3d()

    x_ran = abs(np.diff(x_lim))
    x_mean = np.mean(x_lim)
    y_ran = abs(np.diff(y_lim))
    y_mean = np.mean(y_lim)
    z_ran = abs(np.diff(z_lim))
    z_mean = np.mean(z_lim)

    plot_radius = 0.4 * max([x_ran, y_ran, z_ran])

    # Ajuste para que o grafico não ultrapassar os limites mínimos/maximos
    min_x, max_x = -30, 30
    min_y, max_y = -30, 30
    min_z, max_z = -10, 30

    x0, x1 = x_mean - plot_radius, x_mean + plot_radius
    y0, y1 = y_mean - plot_radius, y_mean + plot_radius
    z0, z1 = z_mean - plot_radius, z_mean + plot_radius

    ax.set_xlim3d([min(x0, min_x), max(x1, max_x)])
    ax.set_ylim3d([min(y0, min_y), max(y1, max_y)])
    ax.set_zlim3d([min(z0, min_z), max(z1, max_z)])


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
