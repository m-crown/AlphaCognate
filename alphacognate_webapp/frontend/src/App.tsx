import { Routes, Route } from "react-router-dom";
import HomePage from "./pages/HomePage";
import AboutPage from "./pages/AboutPage";
import NotFoundPage from "./pages/NotFoundPage";
import { MantineProvider } from '@mantine/core';
import { MobileNavbar } from './components/AppShell';
import StructurePage from "./pages/StructurePage";
import '@mantine/core/styles.css';
import '@mantine/dates/styles.css'; //if using mantine date picker features
import 'mantine-react-table/styles.css'; //make sure MRT styles were imported in your app root (once)

//make a button for this dark theme toggling

function App() {
  return (
    <MantineProvider defaultColorScheme="light">
        <Routes>
        <Route path="/" element={<MobileNavbar />}>
            <Route index element={<HomePage />} />
            <Route path="about" element={<AboutPage />} />
            <Route path="structure/:id" element={<StructurePage />} />
            <Route path="*" element={<NotFoundPage />} />
        </Route>
        </Routes>
      </MantineProvider>
  );
}

export default App;