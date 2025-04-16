import AllStructureTable from "../components/AllStructureTable";
import AutoComplete from "../components/SearchBox";
import { Grid } from "@mantine/core";
import { useNavigate } from "react-router-dom";

function HomePage() {
  const navigate = useNavigate();
  const handleStructureClick = (structure: string) => {
    console.log("Selected Gencode ID:", structure);
    // Navigate to details page with selected row ID as a URL parameter
    navigate(`/structure/${structure}`, { state: { pdbUrl: structure } });
};

  return (
        <Grid>
          <Grid.Col span={12}>
            <AutoComplete />
            <AllStructureTable onStructureRowClick={handleStructureClick}/>
          </Grid.Col>
        </Grid>
  );
}

export default HomePage;
