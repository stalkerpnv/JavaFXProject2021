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

import javafx.scene.control.TextField;
import javafx.stage.Modality;
import javafx.stage.Stage;
import sample.pursuer.Pursuer;
import sample.pursuer.Segment;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.ResourceBundle;

public class SegmentDialogController {
    @FXML
    TextField idS;

    @FXML
    TextField tS;

    @FXML
    TextField psiS;

    private Segment segment;

    private ObservableList<Segment> appMaintvSegments = FXCollections.observableArrayList();
    private ArrayList<Segment> appMainListSegments = new ArrayList<>();

    public void setParameters(int i) {
        idS.setText(Integer.toString(i));
        tS.setText("1000.0");
        psiS.setText("0.0");
    }

    public void addSegment(ActionEvent event) {
        int id = Integer.parseInt(idS.getText().trim());
        double t = Double.parseDouble(tS.getText().trim());
        double psi = Double.parseDouble(psiS.getText().trim());
        segment = new Segment(id, t, psi);
        appMaintvSegments.add(segment);
        appMainListSegments.add(segment);
        closeStage(event);
    }
    public void setObservableListSegment(ObservableList<Segment> tvSegments) {
        this.appMaintvSegments = tvSegments;
    }

    public void setSegmentList(ArrayList<Segment> segmentList) {
        this.appMainListSegments = segmentList;
    }

    private void closeStage(ActionEvent event) {
        Node source = (Node) event.getSource();
        Stage stage = (Stage) source.getScene().getWindow();
        stage.close();
    }
}
