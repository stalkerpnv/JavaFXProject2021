package sample;

import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Node;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.event.ActionEvent;

import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.stage.Modality;
import javafx.stage.Stage;
import sample.pursuer.Pursuer;
import sample.pursuer.Segment;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.ResourceBundle;


public class PursuerDialogController implements Initializable {
    @FXML
    private TextField tfid;
    @FXML
    private TextField tfxc;
    @FXML
    private TextField tfyc;

    @FXML
    private TableView<Segment> tableS;
    @FXML
    private TableColumn colidS;
    @FXML
    private TableColumn colTS;
    @FXML
    private TableColumn colAngleS;
    @FXML
    private Button saveButton;

    private Pursuer pursuer;

    private ObservableList<Pursuer> appMaintvPursuers = FXCollections.observableArrayList();
    private ArrayList<Pursuer> appMainListPursuers = new ArrayList<>();

    private ObservableList<Segment> tvSegments = FXCollections.observableArrayList();
    private ArrayList<Segment> segmentsList = new ArrayList<>();
    private static double l, alpha;

    @FXML
    public void addPursuer(ActionEvent event) {
        int id = Integer.parseInt(tfid.getText().trim());
        double xc = Double.parseDouble(tfxc.getText().trim());
        double yc = Double.parseDouble(tfyc.getText().trim());
        pursuer = new Pursuer(id, xc, yc, l, alpha, segmentsList);
        appMaintvPursuers.add(pursuer);
        appMainListPursuers.add(pursuer);
/*        System.out.println("Add " + pursuer);*/
        closeStage(event);
    }

    public void setObservableList(ObservableList<Pursuer> tvPursuers) {
        this.appMaintvPursuers = tvPursuers;
    }

    public void setPursuerList(ArrayList<Pursuer> pursuersList) {
        this.appMainListPursuers = pursuersList;
    }

    public void setParameters(int index, double xc, double yc, double MainL, double MainAlpha) {
        tfid.setText(Integer.toString(index));
        tfxc.setText(Double.toString(xc));
        tfyc.setText(Double.toString(yc));
        l = MainL;
        alpha = MainAlpha;
    }

    @FXML
    void openDialogSegmentAdd(ActionEvent event) throws IOException {
        FXMLLoader fxmlLoader = new FXMLLoader(getClass().getResource("SegmentDialog.fxml"));
        Parent parent = fxmlLoader.load();
        SegmentDialogController sdc = fxmlLoader.<SegmentDialogController>getController();

        sdc.setParameters(segmentsList.size() + 1);
        sdc.setSegmentList(segmentsList);
        sdc.setObservableListSegment(tvSegments);

        Scene scene = new Scene(parent);
        Stage stage = new Stage();
        stage.setTitle("Add Segment");
        stage.initModality(Modality.APPLICATION_MODAL);
        stage.setScene(scene);
        stage.setResizable(false);
        stage.showAndWait();
        saveButton.setVisible(true);
    }

    private void closeStage(ActionEvent event) {
        Node source = (Node) event.getSource();
        Stage stage = (Stage) source.getScene().getWindow();
        stage.close();
    }

    @Override
    public void initialize(URL location, ResourceBundle resources) {
        colidS.setCellValueFactory(new PropertyValueFactory<>("id"));
        colTS.setCellValueFactory(new PropertyValueFactory<>("t"));
        colAngleS.setCellValueFactory(new PropertyValueFactory<>("angle"));
        tableS.setItems(tvSegments);
        tableS.setPlaceholder(new Label("List of segments is empty "));
        saveButton.setVisible(false);
    }
}

