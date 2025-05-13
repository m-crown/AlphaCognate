
import { useState, useEffect } from 'react';
import { Table, Pagination } from '@mantine/core';

const elements = [
  { position: 1, name: 'Hydrogen', symbol: 'H', mass: 1.008 },
  { position: 2, name: 'Hydrogen', symbol: 'H1', mass: 1.009 },
  { position: 3, name: 'Helium', symbol: 'He', mass: 4.0026 },
  { position: 4, name: 'Carbon', symbol: 'C', mass: 12.011 },
  { position: 5, name: 'Carbon', symbol: 'C1', mass: 12.010 },
  { position: 6, name: 'Carbon', symbol: 'C2', mass: 12.012 },
  { position: 7, name: 'Nitrogen', symbol: 'N', mass: 14.007 },
  { position: 8, name: 'Oxygen', symbol: 'O', mass: 15.999 },
    { position: 9, name: 'Fluorine', symbol: 'F', mass: 18.998 },
  { position: 10, name: 'Neon', symbol: 'Ne', mass: 20.180 },
  { position: 11, name: 'Sodium', symbol: 'Na', mass: 22.990 },
  { position: 12, name: 'Magnesium', symbol: 'Mg', mass: 24.305 },
  { position: 13, name: 'Aluminium', symbol: 'Al', mass: 26.982 },
  { position: 14, name: 'Silicon', symbol: 'Si', mass: 28.085 },
  { position: 15, name: 'Phosphorus', symbol: 'P', mass: 30.974 },
  { position: 16, name: 'Sulfur', symbol: 'S', mass: 32.06 },
  { position: 17, name: 'Chlorine', symbol: 'Cl', mass: 35.45 },
  { position: 18, name: 'Argon', symbol: 'Ar', mass: 39.948 },
  { position: 19, name: 'Potassium', symbol: 'K', mass: 39.098 },
  { position: 20, name: 'Calcium', symbol: 'Ca', mass: 40.078 },
];

const pageSize = 5;

function AboutPage() {
  const [page, setPage] = useState(1);
  const [selectedRow, setSelectedRow] = useState<number | null>(null);

  const paginatedData = elements.slice((page - 1) * pageSize, page * pageSize);

  useEffect(() => {
    const row = elements.find((el) => el.position === selectedRow);
    if (row) console.log('Selected row:', row);
  }, [selectedRow]);

  // Count how many times each name appears
  const nameCounts: Record<string, number> = {};
  for (const el of paginatedData) {
    nameCounts[el.name] = (nameCounts[el.name] || 0) + 1;
  }

  const seenNames = new Set<string>();
  const rows = paginatedData.map((el) => {
    const isFirst = !seenNames.has(el.name);
    if (isFirst) seenNames.add(el.name);

    return (
      <Table.Tr
        key={el.position}
        onClick={() => setSelectedRow(el.position)}
        style={{ cursor: 'pointer' }}
        bg={selectedRow === el.position ? 'var(--mantine-color-blue-light)' : undefined}
      >
        <Table.Td>{el.position}</Table.Td>
        {isFirst && <Table.Td rowSpan={nameCounts[el.name]}>{el.name}</Table.Td>}
        {/* Skip name cell for duplicates */}
        {!isFirst && null}
        <Table.Td>{el.symbol}</Table.Td>
        <Table.Td>{el.mass}</Table.Td>
      </Table.Tr>
    );
  });

  return (
    <>
      <Table>
        <Table.Thead>
          <Table.Tr>
            <Table.Th>Position</Table.Th>
            <Table.Th>Name</Table.Th>
            <Table.Th>Symbol</Table.Th>
            <Table.Th>Mass</Table.Th>
          </Table.Tr>
        </Table.Thead>
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>

      <Pagination
        total={Math.ceil(elements.length / pageSize)}
        value={page}
        onChange={setPage}
        mt="md"
      />
    </>
  );
}

export default AboutPage;

